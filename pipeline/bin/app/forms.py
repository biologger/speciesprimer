#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileRequired
from wtforms import StringField, BooleanField, SubmitField
from wtforms.fields import FieldList, SelectMultipleField, SelectField, IntegerField, FloatField
from wtforms.fields.html5 import EmailField
from wtforms.validators import DataRequired, Email, Optional, Length


class LoginForm(FlaskForm):
    email = EmailField('Email address', validators=[DataRequired(), Email()])
    submit = SubmitField('Sign In')

class PipeOptions(FlaskForm):
    new = SubmitField("New")
    search = SubmitField("Search")
    control = SubmitField("Control")
    change = SubmitField("Change")

class StartForm(FlaskForm):
    search_path = StringField("Path to search for configuration files", default = "/home/primerdesign", validators=[DataRequired()])
    selected_species = FieldList(StringField(""), min_entries=1)
    submit = SubmitField("Search configuration files")

class NewForm(FlaskForm):
    selected_species = FieldList(StringField("", validators=[Length(min=2)]), min_entries=1, validators=[DataRequired()])
    wd_path = StringField("Path to working directory", default = "/home/primerdesign", validators=[DataRequired()])
    same_setting = BooleanField('Use the same settings for all targetspecies')
    submit = SubmitField("Settings")

class ControlRunForm(FlaskForm):
    submit = SubmitField("Start Primerdesign")
    stop = SubmitField("Stop Primerdesign")

class SettingsForm(FlaskForm):
    targets = SelectField(label="Select target", choices=[])
    skip_tree = BooleanField(label="Skip the core gene alignment and tree for visualization and troubleshooting", default = False)
    work_offline = BooleanField("Work offline with local genome assemblies", default = False)
    skip_download = BooleanField("Skip the download of Genomes from NCBI", default = False)
    assemblylevel = SelectMultipleField(
        "Assembly level", choices=[
            ('all', "All"), ('complete', 'Complete'), ('chromosome', 'Chromosome'),
            ('scaffold', 'Scaffold'), ('contig', 'Contig')], validators=[DataRequired()])
    remoteblast = BooleanField(
            "Use remote option of BLAST+?", default = False)
    blastseqs = SelectField("Maximal number of sequences per BLAST search", coerce=int, choices=[(100, "100"), (500, "500"), (1000, "1000"), (2000, "2000"), (5000, "5000")], default=1000)
    qc_genes = SelectMultipleField(
            "Gene(s) for BLAST search in the initial quality control step",
            choices=[('rRNA', "16S rRNA"), ('tuf', "tuf"), ('recA', "recA"), ('dnaK', "dnaK"), ('pheS', "pheS")], validators=[DataRequired()])
    exception = StringField(
        "Primer binding to this non-target species is tolerated", validators=[Optional()])
    minsize = IntegerField("Minimal Amplicon size", default = 70)
    maxsize = IntegerField("Maximal Amplicon size", default = 200)
    designprobe = BooleanField("Pick internal hybridization oligo")
    mfold = FloatField(
            "ΔG threshold for secondary structures in PCR products"
            " at 60°C calculated by mfold", default = -3.0)
    mpprimer = FloatField(
            "ΔG  threshold for 3'-end primer dimer binding", default = -3.5)
    mfeprimer_threshold = SelectField("MFEprimer threshold for nontarget sequence PPC", coerce=int, choices=[(80, "80"), (85, "85"), (90, "90"), (95, "95"), (100, "100")], default=90)
    ignore_qc = BooleanField(
            "Include genomes that did not pass quality control")
    blastdbv5 = BooleanField("BLAST DB Version 5")
    intermediate = BooleanField("Do not delete intermediate files")
    change_wd = StringField("Change path of the working directory", default = "/home/primerdesign")
    submit = SubmitField("Submit settings")
    reset = SubmitField("Reset page")

class ChangeSettingsForm(FlaskForm):
    targets =  FileField(validators=[FileRequired()])
    submit = SubmitField("Submit")

class PipeConfig(FlaskForm):
    down_list = SubmitField("Download")
    up_list = FileField("Upload species list")
    reset_list = SubmitField("Reset to default")
    down_abbrev = SubmitField("Download")
    up_abbrev = FileField("Upload Genus abbreviations file")
    reset_abbrev = SubmitField("Reset to default")
    down_p3sett = SubmitField("Download")
    up_p3sett = FileField("Upload Primer3 settings file")
    reset_p3sett = SubmitField("Reset to default")
    down_noblast = SubmitField("Download")
    up_noblast= FileField("Upload NO_BLAST.gi")
    reset_noblast = SubmitField("Reset to default")
    submit = SubmitField("Upload")

class DownloadDB(FlaskForm):
    update_blastdb = SubmitField("Start download")
    blastdb_firstpart = IntegerField("First part of the BLAST DB", default=0)
    blastdb_parts = IntegerField("Last part of the BLAST DB", default=70)
    delete = BooleanField("Delete archive files after extraction", default=False)
    update_txids = SubmitField("Update")
    get_blastdb = SubmitField("Start download")

class DBForm(FlaskForm):
    submit = SubmitField("Start BLAST db download")
    stop = SubmitField("Stop BLAST db download")
