#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import io
import sys
import json
import codecs
import pandas as pd
from pathlib import Path
from IPython.display import display
from IPython.display import HTML
from IPython.display import clear_output
from traitlets import traitlets

#sys.path.append('speciesprimer/pipeline')

import speciesprimer
# reload for changes in code
import importlib
importlib.reload(speciesprimer)

from basicfunctions import GeneralFunctions as G
from basicfunctions import HelperFunctions as H

from ipywidgets import widgets, interactive, interactive_output, Layout, VBox, HBox, Label


layout = Layout(width="auto")


class SpeciesPrimerConfiguration(object):
    def __init__(self):
        self.pipedir = Path(speciesprimer.__file__).resolve().parent
        self.dictionaries = Path(self.pipedir, "dictionaries")
        self.default_dictionaries = Path(self.dictionaries, "default")
        self.uploaded_settings = {}

    def select_dict_dir(self):
        pass

    def help_for_files(self, key):
        h0 = ("""
        The species list is used to evaluate the
        specificity of the primer. Reported primer
        pairs should not amplify genomic regions of
        species in this list.
        (The evaluation is based on speciesnames in
        the BLAST database)
        """)
        h1 = ("""
        The genus abbreviations file is a .csv file
        with abbreviations for genus names.
        The abbreviations make the names of the
        output files and primer pairs shorter.
        e.g. Lb for Lactobacillus
        """)

        h2 = ("""
        Primer3 settings are saved in this file.
        For an explanation of the parameters see
        https://primer3.ut.ee/primer3web_help.htm
        """)

        h3 = """
        BLAST databases can sometimes contain sequences
        with wrong taxonomic labels. This can prevent
        primer design by speciesprimer because no species
        specific sequences are identified or no specific
        primers can be found. Therefore such sequences can
        be excluded from the specificity evaluation.
        """

        info_dict = {
            'species_list.txt': h0,
             'genus_abbrev.csv': h1,
             'p3parameters': h2,
             'no_blast.gi': h3}

        tips = """
        Click "Reset settings" to display the default settings.
        Click "Save settings" to save the default settings.

        Upload your changed files and click on "Save settings"
        to save the data.

        Click "Display settings" to see the data that will
        be saved.
        """

        print(info_dict[key])
        print(tips)

    def on_selection_change(self, change):
        accept_dict = {
            'species_list.txt': '.txt',
             'genus_abbrev.csv': '.csv',
             'p3parameters': '',
             'no_blast.gi': '.gi'}
        self.upwid.accept = accept_dict[change.new]
        with self.output:
            self.output.clear_output()
            self.help_for_files(change.new)

    def pretty_setting_display(self, key, content):
        if key.endswith(".csv"):
            df = pd.read_csv(io.BytesIO(content), index_col=0, header=None)
            df.index.name = "Genus"
            df.columns = ["Abbreviation"]
            display(df)
        else:
            print(codecs.decode(content, encoding="utf-8"))

    def on_upload_change(self, change):
        with self.output:
            self.output.clear_output()
            key = self.selection.value
            filename = list(change.new.keys())[0]
            content = change.new[filename]["content"]
            self.uploaded_settings.update({key: content})
            print("Uploaded", filename)
            print()
            self.pretty_setting_display(key, content)

    def on_display_request(self, event):
        with self.output:
            self.output.clear_output()
            key = self.selection.value
            if key in self.uploaded_settings:
                content = self.uploaded_settings[key]
                self.pretty_setting_display(key, content)
            else:
                print("Did not find any settings in file upload")

    def save_settings_to_file(self, event):
        with self.output:
            self.output.clear_output()
            dp = self.dictionaries
            key = self.selection.value
            fp = Path(dp, key)
            if key in self.uploaded_settings:
                content = self.uploaded_settings[key]
                with open(fp, "wb") as f:
                    f.write(content)
                print("Saved", key)
            else:
                print("Did not find any changed settings")

    def get_default_settings(self, event):
        with self.output:
            self.output.clear_output()
            dp = self.default_dictionaries
            key = self.selection.value
            fp = Path(dp, key)
            with open(fp, "rb") as f:
                content = f.read()
            self.uploaded_settings.update({key: content})
            self.pretty_setting_display(key, content)

    def pipeline_configuration(self):
        self.output = widgets.Output()
        options = [
            "species_list.txt", "genus_abbrev.csv",
            "p3parameters", "no_blast.gi"]
        self.selection = widgets.Select(
            options=options, value=options[0], rows=len(options))
        self.upwid = widgets.FileUpload(
            multiple=False, accept="", description="Upload file")
        reswid = widgets.Button(
            description='Reset settings',
            style={'description_width': 'initial'})
        displaywid = widgets.Button(
                description='Display settings',
                style={'description_width': 'initial'})
        savewid = widgets.Button(
                description='Save settings',
                style={'description_width': 'initial'})
        label=Label("Alternative path for settings:")
        alt_path = widgets.Text(value="/pipeline/dictionaries")

        r0 = self.selection
        r1 = HBox([self.upwid, savewid])
        r2 = HBox([reswid, displaywid])
        r3 = VBox([label, alt_path])
        controls = VBox([r0, r1, r2, ])#r3])
        applet = HBox([controls, self.output])

        accordion = widgets.Accordion(
            children=[applet],
            titles=["SpeciesPrimer configuration"],
            layout=layout)
        accordion.set_title(0, "SpeciesPrimer configuration")

        display(accordion)

        self.selection.observe(self.on_selection_change, names='value')
        self.upwid.observe(self.on_upload_change, names='value')

        displaywid.on_click(self.on_display_request)
        reswid.on_click(self.get_default_settings)
        savewid.on_click(self.save_settings_to_file)

class WidgetDesigner(object):
    # Create widgets from pipeline/dictionaries/default/batchassist_inputdict.json
    def __init__(self):
        self.pipedir = Path(speciesprimer.__file__).resolve().parent
        self.dictionaries = Path(self.pipedir, "dictionaries")

    def create_widgets_from_batchassist(self):
        text_list = ["customdb", "exception", "path", "strains"]
        widget_dict = {}
        raw_widgets = {}

        config_def = os.path.join(self.dictionaries, "default", "batchassist_inputdict.json")
        with open(config_def) as f:
            for line in f:
                defdict = json.loads(line)

        for key in defdict.keys():
            if "eval" in defdict[key].keys():
                boxed = True
                eval_val = defdict[key]["eval"]
                label = defdict[key]["prompt"].split("\n")[0]

                if eval_val[0] == "boolean" or eval_val[0] == "condboolean":
                    wid = widgets.Checkbox(
                                            description="", #key,
                                            layout=Layout(width='auto',
                                            margin = '0px 0px 5px -80px'))
                    raw_widgets.update({key: wid})
                    boxed = False

                elif eval_val[0] == "number":
                    if eval_val[1][0] == 'int':
                        if "size" in key:
                            slider = widgets.IntSlider(
                                min=50, max=500, value=eval_val[1][1], step=1)
                            text = widgets.BoundedIntText(
                                value=eval_val[1][1], min=50, max=500)
                        else:
                            slider = widgets.IntSlider(value=eval_val[1][1], step=1)
                            text = widgets.IntText(value=eval_val[1][1])

                    else:
                        if eval_val[1][1] < 0:
                            slider = widgets.FloatSlider(
                                value=eval_val[1][1], min=-5, max=0, step=0.1)
                            text = widgets.BoundedFloatText(
                                value=eval_val[1][1], min=-5, max=0)

                        else:
                            slider = widgets.FloatSlider(
                                min=0.000001, max=500,
                                value=eval_val[1][1], step=0.001)
                            text = widgets.BoundedFloatText(
                                value=500.0, max=500.0, min=1e-06)

                    slider.layout = layout
                    widgets.jslink((slider, 'value'), (text, 'value'))
                    raw_widgets.update({key: slider})
                    wid = VBox([slider, text])

                elif eval_val[0] == "options":
                    options=eval_val[1]
                    if "offline" in options:
                        options=eval_val[1][0:-1]
                    wid = widgets.SelectMultiple(
                        options=options, value=[eval_val[1][0]], rows=len(options))
                    raw_widgets.update({key: wid})

                elif eval_val[0] == "option":
                    wid = widgets.Dropdown(
                        options=eval_val[1], value=eval_val[1][0])
                    raw_widgets.update({key: wid})

                else:
                    print("Unknown option")
                    print(eval_val, label)

                if boxed:
                    wid.layout=layout
                    box = VBox([Label(label), wid])
                    box.layout=layout
                    widget_dict.update({key: box})
                else:
                    box = HBox([wid, Label(label)])
                    box.layout=layout
                    widget_dict.update({key: box})

            else:
                if key in text_list:
                    fulllabel = defdict[key]["prompt"]
                    label = fulllabel.split(" or type")[0]
                    label = fulllabel.split(" or hit return")[0]
                    wid = widgets.Textarea()
                    wid.layout = layout
                    raw_widgets.update({key: wid})
                    box = VBox([Label(label), wid])
                    box.layout=layout
                    widget_dict.update({key: box})

        return widget_dict, raw_widgets

    def organize_widgets(self):
        # create widgets
        widget_dict, raw_widgets = self.create_widgets_from_batchassist()

        # add missing options (to do?: include in batchassist options)
        missing_keys = {
            "strains": ["", "Select strains for runmode strain"] ,
            "path":["/primerdesign", "Select path (working directory)"]}

        for k in missing_keys.keys():
            wid = widgets.Textarea(value=missing_keys[k][0])
            wid.layout=layout
            raw_widgets.update({k: wid})
            box = VBox([Label(missing_keys[k][1]), wid])
            widget_dict.update({k: box})

        # prepare widget order and visualization
        readable_settings_dict = {
            "Species specificity settings:":
                ["runmode", "strains", "exception", "virus", ],
            "Download settings:":
                ["assemblylevel", "offline", "skip_download", "genbank"],
            "Input Quality Control settings:":
                ["qc_gene", "ignore_qc"],
            "Computing ressources settings:":
                ["path", "customdb", "blastseqs", "intermediate", "skip_tree", "nolist"],
            "Primer design settings:":
                ["minsize", "maxsize", "probe"],
            "Primer Quality Control settings:":
                ["mpprimer", "mfethreshold", "mfold"],
            "Advanced fine tuning:":
                ["evalue", "nuc_identity"]}

        settings_order = []
        widget_list = []
        for k, vals in readable_settings_dict.items():
            widget_group = []
            for v in vals:
                settings_order.append(v)
                widget_group.append(widget_dict[v])
            display_group = VBox(widget_group)
            widget_list.append(display_group)
        titles = [k for k in readable_settings_dict.keys()]

        return settings_order, raw_widgets, widget_list, titles

class SettingsStorage(object):
    def __init__(
            self,
            settings_order=[], raw_widgets=[],
            widget_list=[], titles=[]):
        self.specieslist = []
        self.settings_dict = {}
        self.settings_order = settings_order
        self.comp_settings_df = pd.DataFrame()
        self.raw_widgets = raw_widgets
        self.widget_list = widget_list
        self.titles = titles
        self.existing_settings = {}
        # get default values from widgets
        self.default_values = {k: v.value for k, v in self.raw_widgets.items()}

    def read_config(self, filepath):
        with open(filepath) as f:
            for line in f:
                setting = json.loads(line)
        return setting

    def write_json_config(self, fp, settings):
        with open(fp, "w") as f:
            f.write(json.dumps(settings))

    def write_config(self, skip_keys=[]):
        d = self.settings_dict
        for k in d.keys():
            if k not in skip_keys:
                fp = self.get_config_path(k)
                print("Write config for\n", k + "\n", fp)
                if not fp.parent.is_dir():
                    os.makedirs(fp.parent)
                self.write_json_config(fp, d[k])
                self.settings_dict.update({k: d[k]})

    def get_settings_from_file(self, config_files):
        self.settings_dict = {}
        for fp in config_files:
            key = fp.parent.parent.name
            config = self.read_config(fp)
            self.settings_dict.update({key: config})

    def get_speciesdir_name(self, k):
        sp_dir = "_".join(k.split(" "))
        if "." in sp_dir:
            sp_dir = "".join(sp_dir.split("."))
        return sp_dir

    def get_config_path(self, key):
        d = self.settings_dict
        sp_dir = self.get_speciesdir_name(key)
        fp = Path(d[key]['path'], sp_dir, "config", "config.json")
        return fp

    def reset_settings(self, event):
        for k, v in self.raw_widgets.items():
            if k in self.default_values.keys():
                v.value = self.default_values[k]

    def load_settings(self, key):
        for k, v in self.raw_widgets.items():
            if k in self.settings_dict[key].keys():
                v.value = self.settings_to_widgets(
                    k, self.settings_dict[key][k])
            elif k in self.default_values.keys():
                v.value = self.default_values[k]
            else:
                v.value = ""

    def settings_to_widgets(self, key, value):
        if key in ["exception", "strains"]:
            value = str(", ".join(value))
        elif isinstance(value, list):
            value = tuple(value)
        else:
            if value is None:
                value = ""
        return value

    def find_old_config_files(self):
        d = self.settings_dict
        for k in d.keys():
            sp_dir = self.get_speciesdir_name(k)
            fp = self.get_config_path(k)
            if fp.is_file():
                print("Config file exists for: " + k)
                print(fp)
                self.existing_settings.update({k: self.read_config(fp)})

    def search_configfile_paths(self, searchpath, filterlist=[]):
        conffiles = []
        limited_dirs = []
        for sp in filterlist:
            limited_dirs.append(self.get_speciesdir_name(sp))
        try:
            dir_path = Path(searchpath)
            for child in dir_path.iterdir():
                if child.is_dir():
                    configfile = Path(child, "config", "config.json")
                    if configfile.is_file():
                        if len(limited_dirs) > 0:
                            if child.name in limited_dirs:
                                conffiles.append(configfile)
                        else:
                            conffiles.append(configfile)
            return conffiles
        except Exception:
            print("No valid path")
            return []

    def settings_dict_to_df(self):
        settings_df = pd.DataFrame(index=self.settings_order)
        d = self.settings_dict
        for k in d.keys():
            settings_df[k] = [
                d[k][v] if v in d[k].keys() else "Not specified"
                for v in self.settings_order]
        return settings_df

    def compare_settings_to_df(self):
        d = self.settings_dict
        comp = [k for k in self.existing_settings.keys() if k in d]
        col = pd.MultiIndex.from_product(
            [comp, ["new", "old"]],
            names=["Setting", "time"])
        self.comp_settings_df = pd.DataFrame(
            index=self.settings_order, columns=col)
        for k in self.existing_settings.keys():
            if k in d.keys():
                new = [d[k][v] for v in self.settings_order]
                old = [
                    self.existing_settings[k][v]
                    if v in self.existing_settings[k].keys()
                    else "Not specified" for v in self.settings_order ]
                self.comp_settings_df[(k, "new")] = new
                self.comp_settings_df[(k, "old")] = old
                print("Config file exists for: " + k)
                print(self.get_config_path(k))

    def parse_widget_values(self, key):
        value = self.raw_widgets[key].value
        if isinstance(value, tuple):
            value = list(value)
        elif key in ["exception", "strains"]:
            spp_list = []
            spp = value.split(",")
            for i, sp in enumerate(spp):
                x = sp.strip()
                if key == "exception":
                    x = x.capitalize()
                if x != "" and x not in spp_list:
                    spp_list.insert(i, x)
            value = spp_list
        else:
            if value == "":
                value = None
        return value


    def get_settings_from_widgets(self, species):
        values = [self.parse_widget_values(k) for k in self.settings_order]
        val_dict = dict(zip(self.settings_order, values))
        sp_dir = self.get_speciesdir_name(species)
        val_dict.update({"target": sp_dir})
        self.settings_dict.update({sp_dir: val_dict})


class CheckConfigFiles(object):
    def __init__(self, output, configstore):
        self.configstore = configstore
        self.existing_settings = configstore.existing_settings
        self.output = output
        self.button_write = widgets.Button(
                        description='Overwrite settings',
                        style={'description_width': 'initial'}
                    )
        self.button_skip = widgets.Button(
                description='Keep settings',
                style={'description_width': 'initial'}
            )

    def write_button_clicked(self, event):
        self.output.clear_output()
        self.configstore.write_config()
        self.button_write.close()
        self.button_skip.close()

    def skip_button_clicked(self, event):
        self.output.clear_output()
        skip_list = list(self.existing_settings.keys())
        print("Skip writing config files for " + "\n".join(skip_list))
        self.configstore.write_config(skip_list)

    def check_existing_configuration(self):
        self.configstore.find_old_config_files()
        with self.output:
            if self.existing_settings != {}:
                self.output.clear_output()
                self.configstore.compare_settings_to_df()

                display(self.configstore.comp_settings_df)
                display(HBox([self.button_skip, self.button_write]))

                self.button_write.on_click(self.write_button_clicked)
                self.button_skip.on_click(self.skip_button_clicked)
            else:
                self.configstore.write_config()



class Settings(object):
    def __init__(self, output, configstore):
        self.configstore = configstore
        self.setting_species = [sp for sp in self.configstore.specieslist]
        self.report_output = widgets.Output()

        self.button_reset = widgets.Button(
                description='Reset settings',
                style={'description_width': 'initial'})
        self.button_settings = widgets.Button(
                description='Submit settings',
                style={'description_width': 'initial'})
        self.speciesnames = widgets.SelectMultiple(
            options=self.configstore.specieslist,
            value=self.configstore.specieslist,
            description='Species:',
            disabled=False
            )
        self.accordion = widgets.Accordion(children=self.configstore.widget_list)
        for i, title in enumerate(self.configstore.titles):
            self.accordion.set_title(i, title)

    def submit_settings(self, event):
        with self.report_output:
            self.report_output.clear_output()
            for species in self.speciesnames.value:
                self.configstore.get_settings_from_widgets(species)
                self.setting_species.remove(species)

            settings_df = self.configstore.settings_dict_to_df()
            display(settings_df)

            if len(self.setting_species) > 0:
                self.speciesnames.options = self.setting_species
                self.speciesnames.value = [self.setting_species[0]]
            else:
                set_w = [self.speciesnames, self.accordion, self.settings_label]
                for w in set_w:
                    w.close()
                self.report_output.clear_output()
                print("Selected settings")
                display(settings_df)
                CCF = CheckConfigFiles(self.report_output, self.configstore)
                CCF.check_existing_configuration()


    def submit_changed_settings(self, event):
        with self.report_output:
            self.report_output.clear_output()
            self.configstore.settings_dict = {}
            for species in self.speciesnames.value:
                self.configstore.get_settings_from_widgets(species)

            settings_df = self.configstore.settings_dict_to_df()
            display(settings_df)

            set_w = [
              self.speciesnames, self.accordion,
              self.settings_label, self.load_setting,
              self.button_change, self.button_reset, self.button_settings,
              self.label]
            for w in set_w:
                w.close()
            self.report_output.clear_output()
            print("Selected settings")
            display(settings_df)
            CCF = CheckConfigFiles(self.report_output, self.configstore)
            CCF.check_existing_configuration()

    def selection_change(self, change):
        self.label.value = ' for: ' + ", ".join(change["new"])

    def get_user_settings(self):
        self.label = widgets.Label(value=' for ' + ", ".join(self.speciesnames.value))
        self.speciesnames.observe(self.selection_change, names='value')
        self.settings_label = HBox([self.button_reset, self.button_settings, self.label])
        self.button_reset.on_click(self.configstore.reset_settings)
        self.button_settings.on_click(self.submit_settings)
        if len(self.setting_species) > 0:
            self.box = VBox([self.speciesnames, self.accordion, self.settings_label, self.report_output])
            display(self.box)

    def request_setting(self, event):
        sp = list(self.speciesnames.value)
        if len(sp) > 1:
            print("Can only load one setting, select one species")
        else:
            sp_dir = self.configstore.get_speciesdir_name(sp[0])
            self.configstore.load_settings(sp_dir)
            print("Loaded settings for " + str(sp[0]))

    def change_settings(self, event):
        with self.report_output:
            self.report_output.clear_output()
            self.settings_label = HBox([self.button_reset, self.button_settings, self.label])
            self.button_reset.on_click(self.configstore.reset_settings)
            self.button_settings.on_click(self.submit_changed_settings)
            self.load_setting.on_click(self.request_setting)

            display(HTML("<h4>Change settings: </h4>"))
            display(self.accordion)
            display(self.settings_label)
            display(self.load_setting)


    def change_user_settings(self):
        self.button_change = widgets.Button(
                description='Change Settings',
                style={'description_width': 'initial'}
            )
        self.load_setting = widgets.Button(
                description='Load Settings',
                style={'description_width': 'initial'}
            )
        self.label = widgets.Label(value=' for ' + ", ".join(self.speciesnames.value))
        self.speciesnames.observe(self.selection_change, names='value')
        self.settings_label = HBox([self.button_change , self.label])
        self.button_change.on_click(self.change_settings)


        self.box = VBox([self.speciesnames, self.settings_label, self.report_output])
        display(self.box)

class TargetSelection(object):
    def __init__(self, configstore):
        self.configstore = configstore
        self.targets = widgets.Textarea(layout=layout)
        self.target_selection = VBox(
            [Label("Please specify target species (comma separated)"),
             self.targets])
        self.example_input = str(
        "Lactobacillus helveticus,\nLactobacillus delbrueckii,"
        "\nLactococcus lactis subsp. lactis""")

    def submit_species(self, event):
        with self.newoutput:
            self.newoutput.clear_output()
            if self.targets.value == "":
                print("Please provide at least one target species")
                print("\nExample:\n" + self.example_input)
            else:
                specieslist = []
                species = self.targets.value.split(",")
                for i, sp in enumerate(species):
                    x = sp.strip()
                    x = x.capitalize()
                    if x != "" and x not in specieslist:
                        specieslist.insert(i, x)
                self.configstore.specieslist = specieslist
                self.configstore.settings_dict = {}
                sel = Settings(self.newoutput, self.configstore)
                sel.get_user_settings()

    def search_old_targets(self, event):
        with self.output:
            self.output.clear_output()
            if self.targets.value == "":
                print("Search for config files")
                config_files = self.configstore.search_configfile_paths(
                                                    self.searchpath_input.value)
            else:
                specieslist = []
                species = self.targets.value.split(",")
                for i, sp in enumerate(species):
                    x = sp.strip()
                    x = x.capitalize()
                    if x != "" and x not in specieslist:
                        specieslist.insert(i, x)

                config_files = self.configstore.search_configfile_paths(
                                        self.searchpath_input.value, specieslist)

            if config_files != []:
                self.configstore.get_settings_from_file(config_files)
                settings_df = self.configstore.settings_dict_to_df()
                display(self.button_change)
                display(settings_df)

            else:
                print(
                    "No config files found in path:",
                    Path(self.searchpath_input.value))

    def new_targets(self):
        self.newoutput = widgets.Output()
        button_send = widgets.Button(
                        description='Submit',
                        style={'description_width': 'initial'}
                    )
        button_send.on_click(self.submit_species)
        vbox_result = widgets.VBox([button_send, self.newoutput])
        dashboard = VBox([self.target_selection, vbox_result])

        accordion = widgets.Accordion(
            children=[dashboard],
            layout=layout)
        accordion.set_title(0, "New run configuration")
        display(accordion)


    def old_targets(self):
        self.output = widgets.Output()
        label = "Directory to search for config files"
        self.searchpath_input = widgets.Text(value="/primerdesign", layout=layout)
        box = VBox([Label(label), self.searchpath_input])
        button_send = widgets.Button(
                        description='Submit',
                        style={'description_width': 'initial'}
                    )
        self.button_change = widgets.Button(
                description='Change Settings',
                style={'description_width': 'initial'}
            )

        comment = Label("Leave empty to search for all files in the path")
        dashboard = VBox([box, self.target_selection, comment, button_send, self.output])

        accordion = widgets.Accordion(
            children=[dashboard],
            layout=layout)
        accordion.set_title(0, "Search existing configuration")
        display(accordion)

        button_send.on_click(self.search_old_targets)
        self.button_change.on_click(self.change_settings)


    def change_settings(self, event):
        with self.output:
            self.output.clear_output()
            specieslist = [
                " ".join(s.split("_")) for s in
                self.configstore.settings_dict.keys()]
            self.configstore.specieslist = specieslist
            sel = Settings(self.output, self.configstore)
            sel.change_user_settings()