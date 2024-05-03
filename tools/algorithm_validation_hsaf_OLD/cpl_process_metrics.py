"""
Class Features

Name:          cpl_process_metrics
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230719'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import inspect

from lib_utils_metrics import BasicMetrics, ExtendedMetrics
from lib_utils_generic import get_dataset_modes, get_dataset_names
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class coupler process metrics
class CplMetrics:

    # method to initialize class
    def __init__(self, dset_interfaces, metrics_seasonal=False,
                 metrics_anomaly=False, metrics_type='extended', swi_option=False):

        self.dset_interfaces = dset_interfaces
        self.metrics_seasonal = metrics_seasonal
        self.metrics_anomaly = metrics_anomaly
        self.swi_option = swi_option

        if metrics_type == 'extended':
            self.metrics_class = ExtendedMetrics
        elif metrics_type == 'basic':
            self.metrics_class = BasicMetrics
        else:
            logging.error(' ===> Dateset metrics "' + metrics_type + '" is not expected by the metrics driver')
            raise NotImplemented('Case not implemented yet')

    # method to setup process metrics
    def setup_metrics(self):

        # set datasets reference
        dset_name_reference = get_dataset_modes(self.dset_interfaces, dset_mode='ref')

        # set datasets names
        dset_name_list = get_dataset_names(dset_name_reference, self.dset_interfaces)

        # set datasets metrics
        dset_metrics_args_common = self.organize_signature(
            dataset_names=dset_name_list, seasonal_metrics=self.metrics_seasonal, swi_option=self.swi_option)
        # get signature
        dset_metrics_signature = self.inspect_signature()
        # filter and adapt signature
        dset_metrics_args_filter = self.filter_signature(dset_metrics_signature, dset_metrics_args_common)
        # initialize class
        dset_metrics_obj = self.metrics_class(**dset_metrics_args_filter)

        return dset_metrics_obj

    # method to filter signature args
    @staticmethod
    def filter_signature(signature_obj, signature_args_common):

        signature_args_filter = {}
        for signature_param in signature_obj.parameters.values():
            param_name = signature_param.name
            if param_name in list(signature_args_common.keys()):
                param_value = signature_args_common[param_name]
            else:
                param_value = signature_param.default

            signature_args_filter[param_name] = param_value

        return signature_args_filter

    # method to organize class signature
    @staticmethod
    def organize_signature(dataset_names=None, seasonal_metrics=False, swi_option=False):
        args_signature = {
            'result_path': None, 'other_name': 'k1',
            'other_name1': 'k1', 'other_name2': 'k2',
            'dataset_names': dataset_names, "seasonal_metrics": seasonal_metrics,
            'swi_option': swi_option,
        }
        return args_signature

    # method to inspect class signature
    def inspect_signature(self):
        return inspect.signature(self.metrics_class)

    # method to check key availability
    @staticmethod
    def __check_key(dset_obj, dset_key='time_start'):
        if dset_key in list(dset_obj.keys()):
            dset_value = dset_obj[dset_key]
        else:
            logging.error(' ===> Dataset key "' + dset_key + '" is not available in the time obj')
            raise RuntimeError('The key is needed by the algorithm. Please check your settings file')
        return dset_value

# -------------------------------------------------------------------------------------
