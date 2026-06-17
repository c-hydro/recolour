#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Automatic launcher for the OpenLayers raster viewer.

Usage:
    python3 launcher.py
    python3 launcher.py --no-browser
    python3 launcher.py --no-venv
    python3 launcher.py --reinstall-deps
    python3 launcher.py --metadata metadata_smap.json
    python3 launcher.py --metadata metadata_smap.json --raster /absolute/path/product.tif

The launcher:
- changes working directory to the folder containing this file;
- creates a local lightweight Python virtual environment when needed;
- installs the dependencies required by render_assets.py;
- regenerates PNG layer overlays and colorbars from the selected metadata file;
- finds a free TCP port automatically;
- starts a local HTTP server;
- opens the default browser at the viewer URL.
"""

from __future__ import annotations

import argparse
import functools
import os
import socket
import subprocess
import sys
import threading
import venv
import webbrowser
from urllib.parse import quote
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path


VIEWER_VERSION = "v12.0-categorical-export-legend"
REQUIRED_IMPORTS = ["numpy", "PIL", "matplotlib", "rasterio"]
VENV_DIRNAME = ".viewer_venv"


class QuietHTTPRequestHandler(SimpleHTTPRequestHandler):
    """Simple handler with compact logs."""

    def log_message(self, format, *args):  # noqa: A003 - keep stdlib signature
        sys.stdout.write("[%s] %s\n" % (self.log_date_time_string(), format % args))


def is_port_free(host: str, port: int) -> bool:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        try:
            sock.bind((host, port))
            return True
        except OSError:
            return False


def find_free_port(host: str, start_port: int, end_port: int) -> int:
    for port in range(start_port, end_port + 1):
        if is_port_free(host, port):
            return port
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.bind((host, 0))
        return int(sock.getsockname()[1])


def venv_python(venv_dir: Path) -> Path:
    if os.name == "nt":
        return venv_dir / "Scripts" / "python.exe"
    return venv_dir / "bin" / "python"


def run(cmd: list[str], cwd: Path | None = None) -> None:
    print("$ " + " ".join(cmd))
    subprocess.run(cmd, cwd=str(cwd) if cwd else None, check=True)


def imports_available(python_exe: Path) -> bool:
    code = "\n".join([f"import {name}" for name in REQUIRED_IMPORTS])
    try:
        subprocess.run(
            [str(python_exe), "-c", code],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
        )
        return True
    except subprocess.CalledProcessError:
        return False
    except FileNotFoundError:
        return False


def ensure_venv(viewer_dir: Path, reinstall: bool = False) -> Path:
    """Create/use a local venv and install only the packages needed for rendering."""
    venv_dir = viewer_dir / VENV_DIRNAME
    py = venv_python(venv_dir)

    if not py.exists():
        print(f"Creating local virtual environment: {venv_dir}")
        venv.EnvBuilder(with_pip=True, clear=False).create(str(venv_dir))

    req = viewer_dir / "requirements.txt"
    if not req.exists():
        req.write_text("numpy\npillow\nmatplotlib\nrasterio\n", encoding="utf-8")

    if reinstall or not imports_available(py):
        print("Installing/updating viewer rendering dependencies ...")
        run([str(py), "-m", "pip", "install", "--disable-pip-version-check", "-r", str(req)], cwd=viewer_dir)
    else:
        print("Rendering dependencies already available in local venv.")

    return py


def main() -> int:
    parser = argparse.ArgumentParser(description="Launch the OpenLayers viewer on a free local port.")
    parser.add_argument("--host", default="127.0.0.1", help="Host/interface to bind. Default: 127.0.0.1")
    parser.add_argument("--start-port", type=int, default=8000, help="First preferred port. Default: 8000")
    parser.add_argument("--end-port", type=int, default=8100, help="Last preferred port. Default: 8100")
    parser.add_argument("--no-browser", action="store_true", help="Do not open the browser automatically")
    parser.add_argument("--index", default="index.html", help="Viewer entry point. Default: index.html")
    parser.add_argument("--metadata", default="metadata.json", help="Metadata configuration file. Default: metadata.json")
    parser.add_argument("--raster", default=None, help="Override metadata product.filename with this GeoTIFF path. Absolute paths are supported.")
    parser.add_argument("--no-regenerate", action="store_true", help="Do not regenerate PNG layers/colorbars before launch")
    parser.add_argument("--no-venv", action="store_true", help="Use the current Python instead of the local .viewer_venv")
    parser.add_argument("--reinstall-deps", action="store_true", help="Force reinstall/update packages in the local venv")
    args = parser.parse_args()

    viewer_dir = Path(__file__).resolve().parent
    print(f"OpenLayers raster viewer {VIEWER_VERSION}")
    print(f"Launcher path: {Path(__file__).resolve()}")
    index_path = viewer_dir / args.index
    try:
        html_text = index_path.read_text(encoding="utf-8") if index_path.exists() else ""
        if "Viewer v12.0" in html_text:
            print("HTML version check: Viewer v12.0 found")
        else:
            print("WARNING: HTML version check: Viewer v12.0 not found")
    except Exception as exc:
        print(f"WARNING: could not read index version: {exc}")

    if not index_path.exists():
        print(f"ERROR: viewer entry point not found: {index_path}", file=sys.stderr)
        return 1

    metadata_path = (viewer_dir / args.metadata).resolve() if not Path(args.metadata).is_absolute() else Path(args.metadata).resolve()
    if not metadata_path.exists():
        print(f"ERROR: metadata file not found: {metadata_path}", file=sys.stderr)
        return 1
    try:
        metadata_path.relative_to(viewer_dir)
    except ValueError:
        print("ERROR: metadata file must be inside the viewer directory so the web app can fetch it.", file=sys.stderr)
        print(f"Given: {metadata_path}", file=sys.stderr)
        return 1
    metadata_url = metadata_path.relative_to(viewer_dir).as_posix()
    print(f"Metadata : {metadata_url}")

    raster_override = None
    if args.raster:
        raster_path = Path(args.raster).expanduser()
        if not raster_path.is_absolute():
            # Relative raster paths are resolved from the current shell directory.
            raster_path = (Path.cwd() / raster_path).resolve()
        else:
            raster_path = raster_path.resolve()
        if not raster_path.exists():
            print(f"ERROR: raster override not found: {raster_path}", file=sys.stderr)
            return 1
        raster_override = raster_path
        print(f"Raster override: {raster_override}")

    os.chdir(viewer_dir)

    if args.no_venv:
        render_python = Path(sys.executable)
        print(f"Using current Python: {render_python}")
    else:
        try:
            render_python = ensure_venv(viewer_dir, reinstall=args.reinstall_deps)
        except subprocess.CalledProcessError as exc:
            print("ERROR: dependency installation failed.", file=sys.stderr)
            print("Check your internet connection or run:", file=sys.stderr)
            print("  python3 -m pip install -r requirements.txt", file=sys.stderr)
            print("Then retry with:", file=sys.stderr)
            print("  python3 launcher.py --no-venv", file=sys.stderr)
            return int(exc.returncode or 1)

    if not args.no_regenerate:
        render_script = viewer_dir / "render_assets.py"
        if render_script.exists():
            print(f"Regenerating layer PNGs and colorbars from {metadata_url} ...")
            try:
                cmd = [str(render_python), str(render_script), "--metadata", metadata_url]
                if raster_override is not None:
                    cmd.extend(["--raster", str(raster_override)])
                run(cmd, cwd=viewer_dir)
            except subprocess.CalledProcessError as exc:
                print(f"ERROR: asset regeneration failed: {exc}", file=sys.stderr)
                return int(exc.returncode or 1)
        else:
            print("WARNING: render_assets.py not found; using existing static PNG files")

    port = find_free_port(args.host, args.start_port, args.end_port)
    url = f"http://{args.host}:{port}/{args.index}?metadata={quote(metadata_url)}"

    handler = functools.partial(QuietHTTPRequestHandler, directory=str(viewer_dir))
    server = ThreadingHTTPServer((args.host, port), handler)

    print("============================================================")
    print(f"OpenLayers viewer launcher ({VIEWER_VERSION})")
    print(f"Directory : {viewer_dir}")
    print(f"Python    : {render_python}")
    print(f"Metadata  : {metadata_url}")
    if raster_override is not None:
        print(f"Raster    : {raster_override}")
    print(f"URL       : {url}")
    print("Stop      : CTRL+C")
    print("============================================================")

    if not args.no_browser:
        threading.Timer(0.7, lambda: webbrowser.open(url)).start()

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping server ...")
    finally:
        server.shutdown()
        server.server_close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
