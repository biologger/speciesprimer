#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileRequired
from wtforms import StringField, BooleanField, SubmitField
from wtforms.fields import FieldList, SelectMultipleField, SelectField
from wtforms.fields import IntegerField, FloatField
from wtforms.fields.html5 import EmailField
from wtforms.validators import DataRequired, Email, Length


class LoginForm(FlaskForm):
    email = EmailField('Email address', validators=[DataRequired(), Email()])
    submit = SubmitField('Sign In')


class PipeOptions(FlaskForm):
    new = SubmitField("New")
    search = SubmitField("Search")
    control = SubmitField("Control")
    change = SubmitField("Change")


class StartForm(FlaskForm):
    search_path = StringField(
        "Path to search for configuration files",
        default="/primerdesign",
        validators=[DataRequired()])
    selected_species = FieldList(StringField(""), min_entries=1)
    submit = SubmitField("Search configuration files")


class NewForm(FlaskForm):
    selected_species = FieldList(StringField(
        "", validators=[Length(min=2)]),
        min_entries=1,
        validators=[DataRequired()])
    wd_path = StringField(
            "Path to working directory",
            default="/primerdesign",
            validators=[DataRequired()])
    same_setting = BooleanField('Use the same settings for all targetspecies')
    submit = SubmitField("Settings")


class ControlRunForm(FlaskForm):
    submit = SubmitField("Start Primerdesign")
    stop = SubmitField("Stop Primerdesign")


class SettingsForm(FlaskForm):
    targets = SelectField(label="Select target", choices=[])
    skip_tree = BooleanField(
        label=(
            "Skip the core gene alignment and tree for visualization "
            "and troubleshooting"),
        default=False)
    offline = BooleanField(
            "Work offline with local genome assemblies", default=False)
    skip_download = BooleanField(
            "Skip the download of Genomes from NCBI", default=False)
    assemblylevel = SelectMultipleField(
        "Assembly level", choices=[
            ('all', "All"),
            ('complete', 'Complete'), ('chromosome', 'Chromosome'),
            ('scaffold', 'Scaffold'), ('contig', 'Contig')],
        default=['all'],
        validators=[DataRequired()])
    customdb = StringField(
            "Specify the path to a custom BLAST database", default=None)

    nolist = BooleanField(
            label=(
                "Species list is not used and only sequences without blast "
                "hits are used for primer design "),
            default=False)

    blastseqs = SelectField(
            "Maximal number of sequences per BLAST search", coerce=int,
            choices=[
                (100, "100"), (500, "500"), (1000, "1000"),
                (2000, "2000"), (5000, "5000")],
            default=1000)
    qc_gene = SelectMultipleField(
            "Gene(s) for BLAST search in the initial quality control step",
            choices=[
                ('rRNA', "16S rRNA"), ('tuf', "tuf"),
                ('recA', "recA"), ('dnaK', "dnaK"), ('pheS', "pheS")],
            default=['rRNA'],
            validators=[DataRequired()])
    exception = FieldList(StringField(""), min_entries=1)
    minsize = IntegerField("Minimal Amplicon size", default=70)
    maxsize = IntegerField("Maximal Amplicon size", default=200)
    probe = BooleanField("Pick internal hybridization oligo", default=False)
    mfold = FloatField(
            "ΔG threshold for secondary structures in PCR products"
            " at 60°C calculated by mfold", default=-3.0)
    mpprimer = FloatField(
            "ΔG  threshold for 3'-end primer dimer binding", default=-3.5)
    mfethreshold = SelectField(
        "MFEprimer threshold for nontarget sequence PPC", coerce=int,
        choices=[
            (80, "80"), (85, "85"), (90, "90"), (95, "95"), (100, "100")],
        default=90)
    ignore_qc = BooleanField(
            "Include genomes that did not pass quality control", default=False)
    intermediate = BooleanField(
            "Do not delete intermediate files", default=False)

    virus = BooleanField(
            "Do you want to design primers for a virus?", default=False)
    genbank = BooleanField(
            "Download genome assemblies from GenBank?", default=False)
    evalue = FloatField(
            "E-value threshold for BLAST search, all results with a lower "
            "value pass.", default=500)
    nuc_identity = IntegerField(
            "Nucleotide identity threshold for BLAST search, all results with "
            "a lower value pass.", default=0)
    runmode = SelectMultipleField(
            "Select runmode", default=['species'],
            choices=[
                ('species', 'species'), 
                ('strain', 'strain'), 
                ('group', 'group')],
            validators=[DataRequired()])
    strains = FieldList(StringField(""), min_entries=1)
    subgroup = FieldList(StringField(""), min_entries=1)
    change_wd = StringField(
            "Change path of the working directory", default="/primerdesign")
    submit = SubmitField("Submit settings")


class ChangeSettingsForm(FlaskForm):
    targets = FileField(validators=[FileRequired()])
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
    up_noblast = FileField("Upload NO_BLAST.gi")
    reset_noblast = SubmitField("Reset to default")
    submit = SubmitField("Upload")


class DownloadDB(FlaskForm):
    update_refprok = SubmitField("Start download ref_prok_rep_genomes")
    update_blastdb = SubmitField("Start download nt")
    delete = BooleanField(
            "Delete archive files after extraction", default=True)
    whichdb = SelectField(
            "Select BLAST DB",
            choices=[
                ("nt", "nt"),
                ("ref_prok_rep_genomes", "ref_prok_rep_genomes")],
            default="ref_prok_rep_genomes")
    get_blastdb = SubmitField("Start download")


class DBForm(FlaskForm):
    submit = SubmitField("Start BLAST db download")
    stop = SubmitField("Stop BLAST db download")
