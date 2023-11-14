"""
Library Features:

Name:          lib_figure_cartography
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230727'
Version:       '1.0.0'
"""


from collections import OrderedDict

import numpy as np

from matplotlib.animation import ArtistAnimation, FFMpegFileWriter
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mc
import matplotlib.patches as mpatches

import cartopy
from cartopy.mpl.geoaxes import GeoAxes
import cartopy.crs as ccrs
from mpl_toolkits.axes_grid1 import AxesGrid

from pygeogrids.grids import BasicGrid, genreg_grid


def plot_map(lons, lats, data, mask=None, title='', title_pad=0,
             suptitle='', cmap=None, vmin=None, vmax=None,
             to_image=False, cmap_n=8, cb_label='', cb_extend='neither',
             llcrnrlat=-65, llcrnrlon=-180,
             urcrnrlat=85, urcrnrlon=180, bg_color='gray',
             cb_xtickloc=None, font_scale=1.4, filename=None,
             cb_xticklabels=None, tmp_dpi=400, figsize=(12, 8),
             fig_dpi=300, fillcontinents=True, cbar_ticks=None,
             grid_loc=(0.05, 0.25, 0.9, 0.65), cbar_pad=0.3):

    data_crs = cartopy.crs.PlateCarree()
    proj = cartopy.crs.PlateCarree()

    axes_class = (GeoAxes, dict(map_projection=proj))

    if isinstance(cmap, str):
        if cmap_n is not None:
            cmap = cm.get_cmap(cmap, cmap_n)
        else:
            cmap = cm.get_cmap(cmap)

    other = BasicGrid(lons, lats)
    reg_grid = genreg_grid(0.1, 0.1, minlat=llcrnrlat, maxlat=urcrnrlat,
                           minlon=llcrnrlon, maxlon=urcrnrlon)

    lut = reg_grid.calc_lut(other, max_dist=25000)
    img = np.ma.masked_where(lut == -1, data[lut])
    cimg = np.flipud(img.reshape(-1, reg_grid.shape[1]))

    fig = plt.figure(figsize=figsize)
    fig.suptitle(suptitle)

    ax = AxesGrid(fig, grid_loc, axes_class=axes_class,
                  nrows_ncols=(1, 1), cbar_location='right',
                  cbar_mode='single', cbar_pad=cbar_pad,
                  cbar_size='3%', label_mode='')

    ax[0].outline_patch.set_visible(False)
    ax[0].set_title(title, y=title_pad)

    # land = cartopy.feature.GSHHSFeature(scale='low', edgecolor='none',
    #                                     facecolor='gray')
    # ax[0].add_feature(land)

    coast = cartopy.feature.NaturalEarthFeature(
        category='physical', scale='50m', name='coastline',
        facecolor='none', edgecolor='k')

    land = cartopy.feature.NaturalEarthFeature(
        category='physical', scale='50m', name='coastline',
        facecolor='gray')

    ax[0].add_feature(land, lw=0.4, zorder=2)
    ax[0].add_feature(coast, lw=0.4, zorder=5)
    ax[0].add_feature(cartopy.feature.BORDERS, lw=0.4, zorder=4)

    img_extent = [llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat]

    norm = mc.Normalize(vmin=vmin, vmax=vmax)
    sc = ax[0].imshow(cimg, norm=norm, cmap=cmap, zorder=3,
                      extent=img_extent, interpolation=None,
                      transform=data_crs)

    cbar = ax.cbar_axes[0].colorbar(sc)
    cbar.set_label_text(cb_label)
    cbar.extend = cb_extend

    if cbar_ticks is not None:
        cbar.cbar_axis.set_ticks(cbar_ticks)

    if filename is not None:
        fig.savefig(filename, dpi=fig_dpi, facecolor='white')


def plot_gmap(lons, lats, data, mask=None, title='', title_pad=0,
              suptitle='', cmap=None, vmin=None, vmax=None,
              to_image=False, cmap_n=8, cb_label='', cb_extend='neither',
              cb_xticks=None, filename=None, norm=None,
              cb_xticklabels=None, tmp_dpi=400, figsize=(12, 6.6), cb_pos=None,
              fig_dpi=300, fillcontinents=True, llcrnrlat=-62, llcrnrlon=-180,
              urcrnrlat=87, urcrnrlon=180, bg_color='gray',
              coastline=None, cbar_size='3%', exp_nc=False, exp_tif=False,
              grid_loc=(0.05, 0.1, 0.9, 0.8), cbar_pad=0.3, show_cbar=True,
              max_dist=25000, proj=None):

    data_crs = cartopy.crs.PlateCarree()

    if proj is None:
        proj = cartopy.crs.PlateCarree()

    # proj = cartopy.crs.Mercator(min_latitude=-65.0)

    axes_class = (GeoAxes, dict(map_projection=proj))

    if isinstance(cmap, str):
        if cmap_n is not None:
            cmap = cm.get_cmap(cmap, cmap_n)
        else:
            cmap = cm.get_cmap(cmap)

    other = BasicGrid(lons, lats)
    reg_grid = genreg_grid(0.1, 0.1, minlat=llcrnrlat, maxlat=urcrnrlat,
                           minlon=llcrnrlon, maxlon=urcrnrlon)

    lut = reg_grid.calc_lut(other, max_dist=max_dist)
    img = np.ma.masked_where(lut == -1, data[lut])
    # cimg = np.flipud(img.reshape(-1, reg_grid.shape[1]))
    cimg = img.reshape(-1, reg_grid.shape[1])

    fig = plt.figure(figsize=figsize)
    fig.suptitle(suptitle)

    if show_cbar:
        ax = AxesGrid(fig, grid_loc, axes_class=axes_class,
                      nrows_ncols=(1, 1), cbar_location='bottom',
                      cbar_mode='single', cbar_pad=cbar_pad,
                      cbar_size=cbar_size, label_mode='')
    else:
        ax = AxesGrid(fig, grid_loc, axes_class=axes_class,
                      nrows_ncols=(1, 1), label_mode='')

    ax[0].outline_patch.set_visible(False)
    ax[0].set_title(title, y=title_pad)

    # land = cartopy.feature.GSHHSFeature(scale='low', edgecolor='none',
    #                                     facecolor='gray')
    # ax[0].add_feature(land)

    coast = cartopy.feature.NaturalEarthFeature(
        category='physical', scale='50m', name='coastline',
        facecolor='none', edgecolor='#777777')

    land = cartopy.feature.NaturalEarthFeature(
        category='physical', scale='50m', name='coastline',
        facecolor=bg_color)

    ax[0].add_feature(land, edgecolor=coastline, lw=0.4, zorder=2)
    # ax[0].add_feature(coast, lw=0.4, zorder=5)
    # ax[0].add_feature(cartopy.feature.BORDERS, lw=0.4, zorder=4)

    img_extent = [llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat]

    if norm is None:
        norm = mc.Normalize(vmin=vmin, vmax=vmax)

    sc = ax[0].imshow(cimg, norm=norm, cmap=cmap, zorder=3,
                      extent=img_extent, interpolation=None,
                      transform=data_crs)

    if show_cbar:

        if cb_pos:
            cax = fig.add_axes(cb_pos)
            ax.cbar_axes[0].set_axis_off()
        else:
            cax = ax.cbar_axes[0]

        cbar = plt.colorbar(sc, cax=cax, extend=cb_extend,
                            orientation='horizontal')
        cbar.set_label(cb_label)

    if cb_xticks is not None:
        cbar.set_ticks(cb_xticks)

    if cb_xticklabels is not None:
        cbar.set_ticklabels(cb_xticklabels)

    if filename is not None:
        fig.savefig(filename, dpi=fig_dpi, facecolor='white')

    if exp_nc:
        import netCDF4
        lat_var = reg_grid.lat2d[:, 0]
        lon_var = reg_grid.lon2d[0, :]

        with netCDF4.Dataset('{:}.nc'.format(filename), 'w') as dst:

            dst.createDimension('lon', lon_var.size)
            dst.createDimension('lat', lat_var.size)

            lon = dst.createVariable('lon', np.float32, ('lon',),
                                     zlib=True, complevel=2)
            lon[:] = lon_var

            lon_attr = OrderedDict({'standard_name': 'longitude',
                                    'long_name': 'location longitude',
                                    'units': 'degrees_east',
                                    'valid_range': (-180.0, 180.0)})
            lon.setncatts(lon_attr)

            lat = dst.createVariable('lat', np.float32, ('lat',),
                                     zlib=True, complevel=2)
            lat[:] = lat_var

            lat_attr = OrderedDict({'standard_name': 'latitude',
                                    'long_name': 'location latitude',
                                    'units': 'degrees_north',
                                    'valid_range': (-90.0, 90.0)})
            lat.setncatts(lat_attr)

            data = dst.createVariable('data', np.float32,
                                      ('lat', 'lon'), zlib=True, complevel=2)

            data[:] = img

    if exp_tif:
        from pyraster.geotiff import GeoTiffFile
        from osgeo import osr

        filename = '{:}.tif'.format(filename)

        lat = reg_grid.lat2d[:, 0]
        lon = reg_grid.lon2d[0, :]

        # set geotransform
        nx = cimg.shape[1]
        ny = cimg.shape[0]

        xmin, ymin, xmax, ymax = lon.min(), lat.min(), lon.max(), lat.max()
        xres = (xmax - xmin) / float(nx)
        yres = (ymax - ymin) / float(ny)
        geotransform = (xmin, xres, 0, ymax, 0, -yres)

        spatialref = osr.SpatialReference()
        spatialref.ImportFromEPSG(4326)
        spatialref_str = spatialref.ExportToWkt()

        geotiff = GeoTiffFile(filename, mode='w', count=1,
                              spatialref=spatialref_str, geotransform=geotransform)
        nodata = -9999.
        cimg.fill_value = nodata

        geotiff.write(np.flipud(cimg.filled()), band=1, nodata=[nodata])
        geotiff.close()

    return sc, ax[0]


def plot_add_gmap(lons, lats, data, ax, cmap=None, cmap_n=None, norm=None, vmin=None, vmax=None,
                  llcrnrlat=-62, llcrnrlon=-180, urcrnrlat=87, urcrnrlon=180, max_dist=15000):

    data_crs = cartopy.crs.PlateCarree()

    if isinstance(cmap, str):
        if cmap_n is not None:
            cmap = cm.get_cmap(cmap, cmap_n)
        else:
            cmap = cm.get_cmap(cmap)

    other = BasicGrid(lons, lats)
    reg_grid = genreg_grid(0.1, 0.1, minlat=llcrnrlat, maxlat=urcrnrlat,
                           minlon=llcrnrlon, maxlon=urcrnrlon)

    if norm is None:
        norm = mc.Normalize(vmin=vmin, vmax=vmax)

    lut = reg_grid.calc_lut(other, max_dist=max_dist)
    img = np.ma.masked_where(lut == -1, data[lut])
    cimg = img.reshape(-1, reg_grid.shape[1])

    img_extent = [llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat]

    ax.imshow(cimg, norm=norm, cmap=cmap, zorder=3,
              extent=img_extent, interpolation=None,
              transform=data_crs)


def plot_gmap_legend(lons, lats, data, mask=None, title='', title_pad=0,
                     suptitle='', cmap=None, vmin=None, vmax=None, to_image=False,
                     cmap_n=8, cb_label='', cb_extend='neither', cb_xticks=None,
                     filename=None, norm=None, cb_xticklabels=None, tmp_dpi=400,
                     figsize=(12, 6.6), fig_dpi=300, fillcontinents=True,
                     llcrnrlat=-65, llcrnrlon=-180, urcrnrlat=85, urcrnrlon=180,
                     bg_color='gray', coastline=None, cbar_size='3%', exp_nc=False,
                     grid_loc=(0.05, 0.1, 0.9, 0.8), cbar_pad=0.3):

    data_crs = cartopy.crs.PlateCarree()
    proj = cartopy.crs.PlateCarree()

    axes_class = (GeoAxes, dict(map_projection=proj))

    if isinstance(cmap, str):
        if cmap_n is not None:
            cmap = cm.get_cmap(cmap, cmap_n)
        else:
            cmap = cm.get_cmap(cmap)

    other = BasicGrid(lons, lats)
    reg_grid = genreg_grid(0.1, 0.1, minlat=llcrnrlat, maxlat=urcrnrlat,
                           minlon=llcrnrlon, maxlon=urcrnrlon)

    lut = reg_grid.calc_lut(other, max_dist=25000)
    img = np.ma.masked_where(lut == -1, data[lut])
    cimg = np.flipud(img.reshape(-1, reg_grid.shape[1]))

    fig = plt.figure(figsize=figsize)
    fig.suptitle(suptitle)

    ax = AxesGrid(fig, grid_loc, axes_class=axes_class,
                  nrows_ncols=(1, 1), label_mode='')

    ax[0].outline_patch.set_visible(False)
    ax[0].set_title(title, y=title_pad)

    # land = cartopy.feature.GSHHSFeature(scale='low', edgecolor='none',
    #                                     facecolor='gray')
    # ax[0].add_feature(land)

    # coast = cartopy.feature.NaturalEarthFeature(
    #     category='physical', scale='50m', name='coastline',
    #     facecolor='none', edgecolor='k')

    land = cartopy.feature.NaturalEarthFeature(
        category='physical', scale='50m', name='coastline',
        facecolor=bg_color)

    ax[0].add_feature(land, edgecolor=coastline, lw=0.4, zorder=2)
    # ax[0].add_feature(coast, lw=0.4, zorder=5)
    # ax[0].add_feature(cartopy.feature.BORDERS, lw=0.4, zorder=4)

    img_extent = [llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat]

    if norm is None:
        norm = mc.Normalize(vmin=vmin, vmax=vmax)

    ax[0].imshow(cimg, norm=norm, cmap=cmap, zorder=3,
                 extent=img_extent, interpolation=None,
                 transform=data_crs)

    patch_list = []
    for label, color in zip(cb_xticklabels, cmap.colors):
        patch_list.append(mpatches.Patch(color=color, label=label))

    ax[0].legend(handles=patch_list, bbox_to_anchor=(
        0., 0, 1., -.102), mode='expand', ncol=9)

    if filename is not None:
        fig.savefig(filename, dpi=fig_dpi, facecolor='white')


def plot_sc(lon, lat, data, vmin=0., vmax=100., cmap='viridis', cb_label=''):

    projection = ccrs.PlateCarree()
    axes_class = (GeoAxes, dict(map_projection=projection))

    fig = plt.figure()
    axgr = AxesGrid(fig, 111, axes_class=axes_class,
                    nrows_ncols=(1, 1), axes_pad=0.6,
                    cbar_location='right', cbar_mode='single',
                    cbar_pad=0.2, cbar_size='3%', label_mode='')

    axgr[0].add_feature(cartopy.feature.LAND)
    axgr[0].add_feature(cartopy.feature.OCEAN)
    axgr[0].add_feature(cartopy.feature.COASTLINE)

    p = axgr[0].scatter(lon, lat, c=data, transform=projection,
                        cmap=cmap, s=2, marker='o', vmin=vmin, vmax=vmax)

    cb = axgr.cbar_axes[0].colorbar(p)
    cb.set_label_text(cb_label)


def plot_im(lon, lat, data, vmin=0., vmax=100., cmap='viridis',
            llcrnrlat=-85, llcrnrlon=-180, urcrnrlat=85, urcrnrlon=180,
            title='', cb_label='', cb_extend='neither', filename=None,
            figsize=(12, 6), cb_size='2%'):

    other = BasicGrid(lon, lat)
    reg_grid = genreg_grid(0.1, 0.1, minlat=llcrnrlat, maxlat=urcrnrlat,
                           minlon=llcrnrlon, maxlon=urcrnrlon)

    lut = reg_grid.calc_lut(other, max_dist=25000)
    img = np.ma.masked_where(lut == -1, data[lut])
    cimg = np.flipud(img.reshape(-1, reg_grid.shape[1]))

    projection = ccrs.PlateCarree()
    axes_class = (GeoAxes, dict(map_projection=projection))

    fig = plt.figure(figsize=figsize)

    grid_loc = (0.05, 0.1, 0.9, 0.8)
    axgr = AxesGrid(fig, grid_loc, axes_class=axes_class,
                    nrows_ncols=(1, 1), axes_pad=0.6,
                    cbar_location='right', cbar_mode='single',
                    cbar_pad=0.18, cbar_size=cb_size, label_mode='')

    axgr[0].set_title(title)

    axgr[0].add_feature(cartopy.feature.LAND)
    axgr[0].add_feature(cartopy.feature.OCEAN)
    axgr[0].add_feature(cartopy.feature.COASTLINE)

    norm = mc.Normalize(vmin=vmin, vmax=vmax)
    img_extent = [llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat]

    p = axgr[0].imshow(cimg, norm=norm, cmap=cmap, zorder=3,
                       extent=img_extent, interpolation=None,
                       transform=projection)

    cb = axgr.cbar_axes[0].colorbar(p, extend=cb_extend)
    cb.set_label_text(cb_label)

    if filename is not None:
        fig.savefig(filename)


class CartoMap():

    def __init__(self, figsize=(12, 6.5), suptitle='', map_crs=None,
                 data_crs=None, map_pos=None, frame=False):

        self.fig = plt.figure(figsize=figsize)
        self.fig.suptitle(suptitle)

        if map_crs is None:
            map_crs = cartopy.crs.PlateCarree()

        if map_pos is None:
            map_pos = [0.05, 0.15, 0.9, 0.8]

        self.ax = self.fig.add_axes(map_pos, projection=map_crs)

        if not frame:
            self.ax.outline_patch.set_visible(False)

        self.coast = cartopy.feature.NaturalEarthFeature(
            category='physical', scale='50m', name='coastline',
            facecolor='none', edgecolor='#222222')

        self.ocean = cartopy.feature.NaturalEarthFeature(
            category='physical', scale='50m', name='ocean',
            facecolor='#b1d1fc')

        self.land = cartopy.feature.NaturalEarthFeature(
            category='physical', scale='50m', name='land',
            facecolor='gray')

        self.border = cartopy.feature.NaturalEarthFeature(
            category='cultural', scale='50m', name='admin_0_countries',
            facecolor='none', edgecolor='#555555')

        self.cax = None
        self.prev_img = None
        self.reg_grid = None

        if data_crs is None:
            self.data_crs = cartopy.crs.PlateCarree()
        else:
            self.data_crs = data_crs

    def plot_anim(self, lon, lat, data, data_extent, **kwargs):

        artists = []

        if 'vmin' not in kwargs:
            vmin = np.nanmin(data)
        else:
            vmin = kwargs['vmin']

        if 'vmax' not in kwargs:
            vmax = np.nanmax(data)
        else:
            vmax = kwargs['vmax']

        if 'cmap' not in kwargs:
            cmap = plt.get_cmap('viridis')
        else:
            cmap = kwargs['cmap']

        img, self.reg_grid = oversample(lon, lat, data, data_extent)

        if self.prev_img is None:
            self.prev_img = np.zeros_like(img)
            fades = [1.]
        else:
            fades = [1 / 3., 2 / 3., 1.]

        for fade in fades:
            f_img = np.zeros(img.shape)
            alpha_channel = np.zeros(img.shape, dtype=np.float32)

            # prev valid but not next
            f_img[~self.prev_img.mask] = self.prev_img[~self.prev_img.mask].data
            alpha_channel[~self.prev_img.mask] = (1. - fade)

            # next valid but not prev
            f_img[~img.mask] = img[~img.mask].data
            alpha_channel[~img.mask] = fade

            # both not valid, nothing to compute
            alpha_channel[img.mask & self.prev_img.mask] = 0

            # both valid
            alpha_channel[~img.mask & ~self.prev_img.mask] = 1.
            f_img[~img.mask & ~self.prev_img.mask] = self.prev_img[
                ~img.mask & ~self.prev_img.mask].data * (1. - fade) + img[
                    ~img.mask & ~self.prev_img.mask].data * fade

            norm = mc.Normalize(vmin=vmin, vmax=vmax, clip=True)(f_img)
            colors = cmap(norm.data)
            colors[..., -1] = alpha_channel

            artists.append(self.plot_img(colors, data_extent, **kwargs))

        self.prev_img = img

        return artists

    def init_moviewriter(self, filename, fps=12, dpi=200):
        self.moviewriter = FFMpegFileWriter(fps)
        self.moviewriter.setup(self.fig, filename, dpi=dpi)

    def finish_moviewriter(self):
        self.moviewriter.finish()

    def plot_anim2(self, lon, lat, data, data_extent, frame=False, **kwargs):

        if 'vmin' not in kwargs:
            vmin = np.nanmin(data)
        else:
            vmin = kwargs['vmin']

        if 'vmax' not in kwargs:
            vmax = np.nanmax(data)
        else:
            vmax = kwargs['vmax']

        if 'cmap' not in kwargs:
            cmap = plt.get_cmap('viridis')
        else:
            cmap = kwargs['cmap']

        if 'cb_extend' not in kwargs:
            cb_extend = 'neither'
        else:
            cb_extend = kwargs['cb_extend']

        if 'cb_label' not in kwargs:
            cb_label = ''
        else:
            cb_label = kwargs['cb_label']

        if 'cb_ticklabels' not in kwargs:
            cb_ticklabels = None
        else:
            cb_ticklabels = kwargs['cb_ticklabels']

        if 'cb_ticks' not in kwargs:
            cb_ticks = None
        else:
            cb_ticks = kwargs['cb_ticks']

        if 'cb_extend' not in kwargs:
            cb_extend = 'neither'
        else:
            cb_extend = kwargs['cb_extend']

        norm = mc.Normalize(vmin=vmin, vmax=vmax, clip=True)
        img, self.reg_grid = oversample(lon, lat, data, data_extent)

        if self.prev_img is None:
            self.prev_img = np.zeros_like(img)
            fades = [1.]
        else:
            fades = [1 / 3., 2 / 3., 1.]

        for fade in fades:
            f_img = np.zeros(img.shape)
            alpha_channel = np.zeros(img.shape, dtype=np.float32)

            # prev valid but not next
            f_img[~self.prev_img.mask] = self.prev_img[~self.prev_img.mask].data
            alpha_channel[~self.prev_img.mask] = (1. - fade)

            # next valid but not prev
            f_img[~img.mask] = img[~img.mask].data
            alpha_channel[~img.mask] = fade

            # both not valid, nothing to compute
            alpha_channel[img.mask & self.prev_img.mask] = 0

            # both valid
            alpha_channel[~img.mask & ~self.prev_img.mask] = 1.
            f_img[~img.mask & ~self.prev_img.mask] = self.prev_img[
                ~img.mask & ~self.prev_img.mask].data * (1. - fade) + img[
                    ~img.mask & ~self.prev_img.mask].data * fade

            colors = cmap(norm(f_img).data)
            colors[..., -1] = alpha_channel

            self.ax.clear()

            artist = self.ax.imshow(colors, zorder=2.2, cmap=cmap,
                                    extent=data_extent, interpolation=None,
                                    transform=self.data_crs, origin='upper')

            if self.cax is not None:
                self.cax.clear()
                self.cax = None

            self.add_cb(artist, extend=cb_extend, label=cb_label,
                        ticklabels=cb_ticklabels, ticks=cb_ticks)

            self.ax.add_feature(self.land, zorder=2.1)
            self.ax.add_feature(self.ocean, zorder=2.3)
            self.ax.add_feature(self.coast, linestyle='-', lw=0.8, zorder=2.4)
            self.ax.add_feature(self.border, linestyle=':', lw=0.6, zorder=2.4)
            if not frame:
                self.ax.outline_patch.set_visible(False)

            self.moviewriter.grab_frame()

        self.prev_img = img

    def plot_lonlat(self, lon, lat, data, data_extent, to_img=True,
                    max_dist=7000., grid_sampling=0.05, **kwargs):
        """
        Plot map.
        """
        if to_img:
            self.img, self.reg_grid = oversample(
                lon, lat, data, data_extent, grid_sampling=grid_sampling,
                max_dist=max_dist)
            artist = self.plot_img(self.img, data_extent, **kwargs)
        else:
            kwargs['data_extent'] = data_extent
            artist = self._plot_lonlat(lon, lat, data, **kwargs)

        return artist

    def _plot_lonlat(self, lon, lat, data, data_extent=None, data_crs=None,
                     map_crs=None, bg_color='gray', title='', cmap=None,
                     cmap_n=8, norm=None, vmin=None, vmax=None,
                     animated=False, alpha=1, zorder=1, show_cb=False,
                     cb_extend='neither', cb_label='', land=False,
                     ocean=False, coast=False, border=False,
                     picker=False, pickradius=0):

        if data_crs is None:
            data_crs = cartopy.crs.PlateCarree()

        if vmin is None:
            vmin = np.nanmin(data)

        if vmax is None:
            vmax = np.nanmax(data)

        if isinstance(cmap, str):
            if cmap_n is not None:
                cmap = plt.get_cmap(cmap, cmap_n)
            else:
                cmap = plt.get_cmap(cmap)

        if norm is None:
            norm = mc.Normalize(vmin=vmin, vmax=vmax)

        artist = self.ax.scatter(lon, lat, c=data, norm=norm, cmap=cmap,
                                 zorder=zorder, transform=data_crs,
                                 vmin=vmin, vmax=vmax)

        if land:
            self.ax.add_feature(self.land, zorder=2.1)
        if ocean:
            self.ax.add_feature(self.ocean, zorder=2.3)
        if coast:
            self.ax.add_feature(self.coast, linestyle='-', lw=0.8, zorder=2.4)
        if border:
            self.ax.add_feature(self.border, linestyle=':', lw=0.6, zorder=2.4)

        self.ax.set_title(title)

        if data_extent:
            self.ax.set_extent(data_extent)

        if show_cb:
            self.add_cb(artist, extend=cb_extend, label=cb_label)

    def plot_img(self, data, data_extent, data_crs=None, map_crs=None,
                 bg_color='gray', title='', cmap=None, cmap_n=8,
                 norm=None, vmin=None, vmax=None, animated=False, alpha=1,
                 zorder=2.2, show_cb=False, cb_extend='neither', cb_label='',
                 cb_pos=None, cb_labelpad=None, cb_ticklabels=None, cb_ticks=None,
                 land=False, ocean=False, coast=False, border=False,
                 picker=False, pickradius=0,
                 title_fontsize=20, title_fontweight='bold', cb_fontsize=10, cb_fontweight='bold'):
        """
        Plot 2d image.

        Parameters
        ----------
        data : numpy.ndarray
        """
        if data_crs is None:
            data_crs = cartopy.crs.PlateCarree()

        if vmin is None:
            vmin = np.nanmin(data)

        if vmax is None:
            vmax = np.nanmax(data)

        if isinstance(cmap, str):
            if cmap_n is not None:
                cmap = plt.get_cmap(cmap, cmap_n)
            else:
                cmap = plt.get_cmap(cmap)

        if norm is None:
            norm = mc.Normalize(vmin=vmin, vmax=vmax)

        artist = self.ax.imshow(data, norm=norm, cmap=cmap, zorder=zorder,
                                extent=data_extent, interpolation=None,
                                transform=data_crs, origin='upper',
                                vmin=vmin, vmax=vmax, picker=picker)

        if land:
            self.ax.add_feature(self.land, zorder=2.1)
        if ocean:
            self.ax.add_feature(self.ocean, zorder=2.3)
        if coast:
            self.ax.add_feature(self.coast, linestyle='-', lw=0.8, zorder=2.4)
        if border:
            self.ax.add_feature(self.border, linestyle=':', lw=0.6, zorder=2.4)

        # set title
        if title_fontweight is not None:
            self.ax.set_title(title, fontsize=title_fontsize, fontweight=title_fontweight)
        else:
            self.ax.set_title(title, fontsize=title_fontsize)

        if show_cb:
            self.add_cb(artist, pos=cb_pos, extend=cb_extend, label=cb_label,
                        labelpad=cb_labelpad, ticklabels=cb_ticklabels,
                        ticks=cb_ticks, fontsize=cb_fontsize, fontweight=cb_fontweight)

        return artist

    def add_cb(self, p, label=None, pos=None, extend='neither', labelpad=None,
               ticks=None, ticklabels=None, animated=False, fontsize=10, fontweight=None):
        """
        Add colorbar.
        """
        if pos is None:
            pos = [0.3, 0.11, 0.4, 0.02]

        if self.cax is None:
            self.cax = self.fig.add_axes(pos)

        cb = plt.colorbar(p, cax=self.cax, extend=extend, orientation='horizontal')

        if label is not None:
            if fontweight is not None:
                cb.set_label(label=label, labelpad=labelpad, size=fontsize, weight=fontweight)
            else:
                cb.set_label(label=label, labelpad=labelpad, size=fontsize)

        if ticks is not None:
            cb.set_ticks(ticks)

        if ticklabels is not None:
            cb.set_ticklabels(ticklabels)
            cb.ax.tick_params(labelsize=fontsize)

    def export_fig(self, filename, dpi=300):
        """
        Export figure to png/pdf.

        Parameter
        ---------
        filename : str
            Filename of graphic.
        dpi : int
            Resolution.
        """
        self.fig.savefig(filename, dpi=dpi)

    def export_animation(self, filename, artists, fps=2, interval=50):
        """
        Export animation.

        Parameters
        ----------
        filename : str
            Filename of animation.
        artists : list
            List of artists.
        fps : int
            Frames per second.
        interval : int
            Interval.
        """
        animation = ArtistAnimation(
            self.fig, artists, interval=interval, repeat_delay=1000)

        animation.save(filename, writer='imagemagick', fps=fps)

    def export_netcdf(self, filename):
        """
        Export image to NetCDF4.

        Parameters
        ----------
        filename : str
            Filename of netCDF file.
        """
        import netCDF4

        lat_var = self.reg_grid.lat2d[:, 0]
        lon_var = self.reg_grid.lon2d[0, :]

        with netCDF4.Dataset(filename, 'w') as dst:

            dst.createDimension('lon', lon_var.size)
            dst.createDimension('lat', lat_var.size)

            lon = dst.createVariable('lon', np.float32, ('lon',),
                                     zlib=True, complevel=2)
            lon[:] = lon_var

            lon_attr = OrderedDict({'standard_name': 'longitude',
                                    'long_name': 'location longitude',
                                    'units': 'degrees_east',
                                    'valid_range': (-180.0, 180.0)})
            lon.setncatts(lon_attr)

            lat = dst.createVariable('lat', np.float32, ('lat',),
                                     zlib=True, complevel=2)
            lat[:] = lat_var

            lat_attr = OrderedDict({'standard_name': 'latitude',
                                    'long_name': 'location latitude',
                                    'units': 'degrees_north',
                                    'valid_range': (-90.0, 90.0)})
            lat.setncatts(lat_attr)

            data = dst.createVariable('data', np.float32,
                                      ('lat', 'lon'), zlib=True, complevel=2)

            data[:] = self.img

    def export_geotiff(self, filename):
        """
        Export to GeoTiff.

        Parameters
        ----------
        filename : str
            Filename of Geotiff file.
        """
        from pyraster.geotiff import GeoTiffFile
        from osgeo import osr

        nodata = -9999.

        lat = self.reg_grid.lat2d[:, 0]
        lon = self.reg_grid.lon2d[0, :]

        ny, nx = self.img.shape

        xmin, ymin, xmax, ymax = lon.min(), lat.min(), lon.max(), lat.max()
        xres = (xmax - xmin) / float(nx)
        yres = (ymax - ymin) / float(ny)
        geotransform = (xmin, xres, 0, ymax, 0, -yres)

        spatialref = osr.SpatialReference()
        spatialref.ImportFromEPSG(4326)
        spatialref_str = spatialref.ExportToWkt()

        self.img.fill_value = nodata

        geotiff = GeoTiffFile(
            filename, mode='w', count=1, spatialref=spatialref_str,
            geotransform=geotransform)

        geotiff.write(np.flipud(self.img.filled()), band=1, nodata=[nodata])
        geotiff.close()


# method to oversampling the data
def oversample(lon, lat, data, extent, grid_sampling=0.05, max_dist=7000.):
    """
    Oversample data.
    """
    other = BasicGrid(lon, lat)
    reg_grid = genreg_grid(grid_sampling, grid_sampling,
                           minlat=extent[2], maxlat=extent[3],
                           minlon=extent[0], maxlon=extent[1])

    lut = reg_grid.calc_lut(other, max_dist=max_dist)
    img = np.ma.masked_where(lut == -1, data[lut])
    img[np.isnan(img)] = np.ma.masked

    return img.reshape(-1, reg_grid.shape[1]), reg_grid
