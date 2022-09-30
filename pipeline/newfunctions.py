#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 11:23:03 2022

@author: ags-bioinfo
"""

class DirectoryStorage(object):
    def __init__(self):
        appname = "app"
        self.HOMEDIR = Path(
                    str(Path(".").resolve()).split(appname)[0])
        self.CONFIGDIR = Path(self.HOMEDIR, appname, "config")

        yaml_file = Path(self.CONFIGDIR, "app_directories.yml")
        self.set_paths_from_yaml(yaml_file, self, self.HOMEDIR)

    def load_dict_from_yaml(self, yaml_fp):
        with open(yaml_fp) as f:
            data = yaml.load(f, Loader=SafeLoader)
        return data

    def set_paths_from_yaml(self, yaml_fp, object, parentdir):
        data = self.load_dict_from_yaml(yaml_fp)
        dirname_dict = {}
        for pair in self.nested_dict_pairs_iterator(data):
            p = list(pair)
            while "subdirs" in p:
                p.remove('subdirs')
            k = p.index("name")
            path, name, attr_name = p[:k], p[-1], p[(k-1)]
            dirname_dict.update({attr_name: name})
            x = [dirname_dict[s] for s in path]
            attr_name, dirpath = path[-1], Path(parentdir, "/".join(x))
            setattr(object, attr_name, dirpath)

    def nested_dict_pairs_iterator(self, dict_obj):
        for key, value in dict_obj.items():
            if isinstance(value, dict):
                for pair in self.nested_dict_pairs_iterator(value):
                    yield (key, *pair)
            else:
                yield (key, value)
                
def check_valid_kws(kws_dict, function, ignore=False):
    subfunc_dict = {
        'facet_kws': sns.FacetGrid,
        'subplot_kws': plt.subplots,
        'gridspec_kws': plt.GridSpec
        }
    if kws_dict == '':
        kws_dict = {}
    else:
        if not isinstance(kws_dict, dict):
            kws_dict = ast.literal_eval(kws_dict)
        if ignore:
            return kws_dict
        sig = inspect.signature(function)
        valid = {k:v for k,v in kws_dict.items() if k in sig.parameters.keys()}
        kws_dict = valid
        for k in valid.keys():
            if "_kw" in k:
                if k in subfunc_dict.keys():
                    kws_dict[k] = check_valid_kws(valid[k], subfunc_dict[k])
    return kws_dict

def check_valid_kws(self, kws_dict, function, ignore=False):
    subfunc_dict = {
        'facet_kws': sns.FacetGrid,
        'subplot_kws': plt.subplots,
        'gridspec_kws': plt.GridSpec
        }
    if kws_dict == '':
        kws_dict = {}
    else:
        if not isinstance(kws_dict, dict):
            kws_dict = ast.literal_eval(kws_dict)
        if ignore:
            return kws_dict
        sig = inspect.signature(function)
        valid = {k:v for k,v in kws_dict.items() if k in sig.parameters.keys()}
        not_valid = [k for k in kws_dict.keys() if not k in sig.parameters.keys()]
        kws_dict = valid
        for k in valid.keys():
            if "_kw" in k:
                if k in subfunc_dict.keys():
                    kws_dict[k] = self.check_valid_kws(valid[k], subfunc_dict[k])

            if self.debug:
                with self.errorout:
                    if valid == {}:
                        print("No valid " + function.__name__ + " Keywords found")
                    if not_valid != []:
                        print("Ignored keywords for " + function.__name__ + ":\n", not_valid)
        return kws_dict