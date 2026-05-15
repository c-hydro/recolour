#!/bin/bash

ENV_FOLDER="/root/envs/env_validation_tools"

source "${ENV_FOLDER}/bin/activate"

echo "Validation HSAF environment loaded"
echo "Python: $(which python)"
echo "Pip:    $(which pip)"
