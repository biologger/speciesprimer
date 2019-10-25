#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from flask import render_template, flash, redirect, url_for, send_file
from app import app
from app.forms import LoginForm, StartForm, NewForm, PipeConfig, PipeOptions, DownloadDB
from app.forms import SettingsForm, ChangeSettingsForm, ControlRunForm, DBForm
import os
import time
import json
import subprocess
import shutil
import signal
from werkzeug.utils import secure_filename


app_dir = os.path.abspath(__file__)
pipe_dir = app_dir.split('bin')[0]
tmp_db_path = os.path.join(pipe_dir, 'tmp_config.json')

def update_tmp_db(update_dict):
    with open(tmp_db_path, 'r') as f:
        for line in f:
            tmp_db = json.loads(line)
    tmp_db.update(update_dict)
    with open(tmp_db_path, 'w') as f:
        f.write(json.dumps(tmp_db))

@app.route('/', methods=['GET', 'POST'])
def login():
    form = LoginForm()
    if os.path.isfile(tmp_db_path):
        with open(tmp_db_path) as f:
            for line in f:
                tmp_db = json.loads(line)
        try:
            mail = tmp_db['email']
            if '@' and '.' in mail:
                return redirect(url_for('index'))
        except KeyError:
            pass
    if form.validate_on_submit():
        email = form.email.data
        if os.path.isfile(tmp_db_path):
            with open(tmp_db_path) as f:
                for line in f:
                    tmp_db = json.loads(line)
            tmp_db.update({'email': email})
        else:
            tmp_db = {
                'email': email, 'new_run':{
                    'path': '', 'targets':{},
                    'same_settings': False,
                    'change_settings': False}}
        with open(tmp_db_path, 'w') as f:
            f.write(json.dumps(tmp_db))
        flash('Saved E-mail address {}'.format(email))
        return redirect(url_for('index'))
    return render_template('login.html', title='Sign In', form=form)

@app.route('/login', methods=['GET', 'POST'])
def change_login():
    form = LoginForm()
    if form.validate_on_submit():
        email = form.email.data
        if os.path.isfile(tmp_db_path):
            with open(tmp_db_path) as f:
                for line in f:
                    tmp_db = json.loads(line)
            tmp_db.update({'email': email})
        else:
            tmp_db = {'email': email}
        with open(tmp_db_path, 'w') as f:
            f.write(json.dumps(tmp_db))
        flash('Saved E-mail address {}'.format(email))
        return redirect(url_for('index'))
    return render_template('change_login.html', title='Sign In', form=form)

@app.route('/index', methods=['GET', 'POST'])
def index():
    return render_template('index.html', title='Home')

@app.route('/new', methods=['GET', 'POST'])
def new():
    form = NewForm()
    if form.validate_on_submit():
        path = form.wd_path.data
        target_species = form.selected_species.data[0].split(',')
        form_dict = {
            'path': path, 'targets':{}, 'same_settings': form.same_setting.data,
            'change_settings': False, 'modus': 'new'}
        for target in list(target_species):
            target = '_'.join(target.strip(' ').split(' '))
            target = target.capitalize()
            form_dict['targets'].update({target:{}})
        update_tmp_db({'new_run': form_dict})
        return redirect(url_for('settings'))
    return render_template('new.html', title='Start Primerdesign', form=form)

@app.route('/continue', methods=['GET', 'POST'])
def search_continue():
    form = StartForm()
    if form.validate_on_submit():
        path = form.search_path.data
        if form.selected_species.data == ['']:
            target_list = None
        else:
            target_list = []
            target_species = form.selected_species.data[0].split(',')
            for target in target_species:
                target = target.capitalize()
                target_list.append(target)
        update_tmp_db({'new_run':{"targets": target_list, "path": path, 'modus': 'continue'}})
        return redirect(url_for('controlrun'))
    return render_template('continue.html', title='Continue Primerdesign', form=form)

@app.route('/change_settings', methods=['GET', 'POST'])
def change_settings():
    form = ChangeSettingsForm()
    if form.validate_on_submit():
        f = form.targets.data
        filename = secure_filename(f.filename)
        if f and allowed_file(filename):
            for line in f:
                old_settings = json.loads(line.decode("utf-8"))
        path = old_settings['path']
        update_dict =  {"new_run":{
            'path': path, 'same_settings':False,
            'change_settings': True, 'modus': 'new',
            'targets': {old_settings['target']: old_settings}}}
        update_tmp_db(update_dict)
        return redirect(url_for('settings'))
    return render_template('change_settings.html', title='Change settings', form=form)

def check_targets():
    target_list = []
    target_choice_list = []
    with open(tmp_db_path, 'r') as f:
        for line in f:
            tmp_db = json.loads(line)
    if tmp_db['new_run']['change_settings']:
        for key in tmp_db['new_run']['targets']:
            key = key.capitalize()
            target_list.append(key)
    else:
        for key in tmp_db['new_run']['targets']:
            if tmp_db['new_run']['targets'][key] == {}:
                key = key.capitalize()
                target_list.append(key)
        target_list.sort()

    if tmp_db['new_run']["same_settings"]:
        for item in target_list:
            item = item.capitalize()
            target_choice_list.append( ' '.join(item.split('_')))
        return target_choice_list, tmp_db
    else:
        for item in target_list:
            item = item.capitalize()
            target_choice_list.append((item, ' '.join(item.split('_'))))
        return target_choice_list, tmp_db

def update_db(sub_dict, target, data):
    (
        qc_genes, ignore_qc, skip_download, assemblylevel, skip_tree, exception,
        minsize, maxsize, designprobe, mfold, mpprimer, mfeprimer_threshold,
        offline, customdb, blastseqs, path, blastdbv5, intermediate, nolist
    ) = data
    sub_dict.update({
            'target': target,
            'minsize': minsize, 'maxsize': maxsize, 'mpprimer': mpprimer,
            'exception': exception, 'path': path, 'intermediate': intermediate,
            'qc_gene': qc_genes, 'mfold': mfold, 'skip_download': skip_download,
            'assemblylevel': assemblylevel, 'skip_tree': skip_tree,
            'nolist': nolist, 'offline': offline, 'ignore_qc': ignore_qc,
            'mfethreshold': mfeprimer_threshold, 'customdb': customdb,
            'blastseqs': blastseqs, 'probe': designprobe, 'blastdbv5': blastdbv5})
    return sub_dict

def reset_settings(form):
    (
    form.qc_genes.data, form.ignore_qc.data, form.skip_download.data,
    form.assemblylevel.data, form.skip_tree.data, form.exception.data,
    form.minsize.data, form.maxsize.data, form.designprobe.data,
    form.mfold.data,form.mpprimer.data, form.mfeprimer_threshold.data,
    form.work_offline.data, form.customdb.data, form.blastseqs.data,
    form.blastdbv5.data, form.intermediate.data, form.nolist.data
    ) = (
    ["rRNA"], False, False,
    ["all"], False, None,
    70, 200, False,
    -3.0, -3.5, 90,
    False, None, 1000,
    False, False, False)
    return form


def load_settings(tmp_db):
    for key in tmp_db['new_run']['targets']:
        target = key
    settings = tmp_db['new_run']['targets'][target]
    (
    qc_genes, ignore_qc,skip_download,
    assemblylevel, skip_tree, exception,
    minsize, maxsize, designprobe,
    mfold,mpprimer, mfeprimer_threshold,
    work_offline, customdb, blastseqs,
    change_wd, blastdbv5, intermediate, nolist
    ) = (
    settings["qc_gene"], settings["ignore_qc"], settings["skip_download"],
    settings["assemblylevel"], settings["skip_tree"], settings["exception"],
    settings["minsize"], settings["maxsize"], settings["probe"],
    settings["mfold"], settings["mpprimer"], settings["mfethreshold"],
    settings["offline"], settings["customdb"], settings["blastseqs"],
    settings['path'], settings['blastdbv5'], settings["intermediate"],
    settings['nolist']
    )

    data = {
        "targets": target, "qc_genes": qc_genes, "ignore_qc": ignore_qc, "skip_download": skip_download,
        "assemblylevel": assemblylevel, "skip_tree": skip_tree, "exception": exception,
        "minsize": minsize, "maxsize": maxsize, "designprobe": designprobe,
        "mfold": mfold, "mpprimer": mpprimer, "mfeprimer_threshold": mfeprimer_threshold,
        "work_offline": work_offline, "customdb": customdb, "blastseqs": blastseqs,
        "change_wd": change_wd, "blastdbv5": blastdbv5,
        "intermediate": intermediate, "nolist": nolist}
    return data

def get_settings(form):
    (
        qc_genes, ignore_qc, skip_download, assemblylevel, skip_tree, exception,
        minsize, maxsize, designprobe, mfold, mpprimer, mfeprimer_threshold,
        offline, customdb, blastseqs, path, blastdbv5, intermediate, nolist
    ) = (
        form.qc_genes.data, form.ignore_qc.data,form.skip_download.data,
        form.assemblylevel.data, form.skip_tree.data, form.exception.data,
        form.minsize.data, form.maxsize.data, form.designprobe.data,
        form.mfold.data,form.mpprimer.data, form.mfeprimer_threshold.data,
        form.work_offline.data, form.customdb.data, form.blastseqs.data,
        form.change_wd.data, form.blastdbv5.data, form.intermediate.data,
        form.nolist.data
        )
    if offline:
        skip_download = True
        assemblylevel = ['offline']
    if skip_download:
        assemblylevel = ['offline']
    if exception == '':
        exception = None
    if customdb == "":
        customdb = None
    elif exception is not None:
        exception = '_'.join(exception.strip(' ').split(' '))
    if 'all' in assemblylevel:
        assemblylevel = ['all']
    return (
        qc_genes, ignore_qc, skip_download, assemblylevel, skip_tree, exception,
        minsize, maxsize, designprobe, mfold, mpprimer, mfeprimer_threshold,
        offline, customdb, blastseqs, path, blastdbv5, intermediate, nolist)


@app.route('/settings', methods=['GET', 'POST'])
def settings():
    target_choice_list, tmp_db = check_targets()
    if len(target_choice_list) == 0:
        return redirect(url_for('controlrun'))

    if tmp_db['new_run']['same_settings']:
        def_path = tmp_db['new_run']['path']
        form = SettingsForm(change_wd = def_path)
        form.targets.choices = [("All targets", "All targets")]
    elif tmp_db['new_run']['change_settings']:
        form = SettingsForm(data=load_settings(tmp_db))
        form.targets.choices = target_choice_list
    else:
        def_path = tmp_db['new_run']['path']
        form = SettingsForm(change_wd = def_path)
        form.targets.choices = target_choice_list

    if form.validate_on_submit():
        if form.reset.data is True:
            reset_form = reset_settings(form)
            return render_template('settings.html', title='Settings for primerdesign', form=reset_form)

        settings_data = get_settings(form)
        if tmp_db['new_run']['same_settings']:
            for target in tmp_db['new_run']['targets']: ## for same settings
                config = tmp_db['new_run']['targets'][target]
                new_config = update_db(config, target, settings_data)
                tmp_db['new_run']['targets'][target].update(new_config)
                with open(tmp_db_path, 'w') as f:
                    f.write(json.dumps(tmp_db))
        else:
            for target in [form.targets.data]:
                config = tmp_db['new_run']['targets'][target]
                new_config = update_db(config, target, settings_data)
                tmp_db['new_run']['targets'][target].update(new_config)
                with open(tmp_db_path, 'w') as f:
                    f.write(json.dumps(tmp_db))

        if tmp_db['new_run']['same_settings']:
            flash('Saved settings for {}'.format(target_choice_list))
            return redirect(url_for('controlrun'))
        elif tmp_db['new_run']['change_settings']:
            conf_path = tmp_db['new_run']['targets'][target]['path']
            save_to = os.path.join(conf_path, target, "config", "config.json")
            try:
                with open(save_to, 'w') as f:
                    f.write(json.dumps(new_config))
                flash('Changed settings for {}'.format(' '.join(target.split('_'))))
                flash(json.dumps(new_config))
                return redirect(url_for('controlrun'))
            except FileNotFoundError:
                flash("The selected directory does not exist")
                return redirect(url_for('settings'))

        else:
            flash('Saved settings for {}'.format(' '.join(target.split('_'))))
        flash(json.dumps(new_config, sort_keys=True))
        return redirect(url_for('settings'))

    return render_template('settings.html', title='Settings for primerdesign', form=form)

@app.route('/primerdesign', methods=['GET', 'POST'])
def primerdesign():
    form = PipeOptions()
    if form.validate_on_submit():
        if form.new.data:
            return redirect(url_for('new'))
        elif form.search.data:
            return redirect(url_for('search_continue'))
        elif form.control.data:
            return redirect(url_for('controlrun'))
        elif form.change.data:
            return redirect(url_for('change_settings'))

    return render_template('speciesprimer.html', title='Primerdesign', form=form)

def start_pipeline():
    subprocess.Popen(["speciesprimerdaemon.py", "start", "100"])
    today = time.strftime("%Y_%m_%d", time.localtime())
    log_file = os.path.join("/", "home", "primerdesign", "speciesprimer_" + today + ".log")
    if os.path.isfile("/tmp/frontail.pid"):
        with open("/tmp/frontail.pid") as f:
            for line in f:
                pid = line.strip()
        pidint = int(pid)
        if isinstance(pidint, int):
            try:
                os.kill(pidint, signal.SIGTERM)
            except ProcessLookupError:
                pass
        os.remove("/tmp/frontail.pid")
    frontail_cmd = [
            "frontail-linux", "-d", "-n", "20", "--pid-path", "/tmp/frontail.pid", log_file]
    subprocess.Popen(frontail_cmd)

def stop_pipeline():
    subprocess.Popen(["speciesprimerdaemon.py", "stop", "100"])

def start_db_download():
    subprocess.Popen(["getblastdbdaemon.py", "start", "200"])
    today = time.strftime("%Y_%m_%d", time.localtime())
    log_file = os.path.join("/", "home", "primerdesign", "speciesprimer_" + today + ".log")
    if os.path.isfile("/tmp/frontail.pid"):
        with open("/tmp/frontail.pid") as f:
            for line in f:
                pid = line.strip()
        pidint = int(pid)
        if isinstance(pidint, int):
            try:
                os.kill(pidint, signal.SIGTERM)
            except ProcessLookupError:
                pass
        os.remove("/tmp/frontail.pid")
    frontail_cmd = [
            "frontail-linux", "-d", "-n", "20", "--pid-path", "/tmp/frontail.pid", log_file]
    subprocess.Popen(frontail_cmd)

def stop_db_download():
    subprocess.Popen(["getblastdbdaemon.py", "stop", "200"])

@app.route('/dbdownload', methods=['GET', 'POST'])
def dbdownload():
    form = DBForm()
    if form.validate_on_submit():
        if form.submit.data:
            start_db_download()
            time.sleep(1)
        elif form.stop.data:
            stop_db_download()
    return render_template('dbdownload.html', title='Control BLAST db download', form=form)

def start_updatedb():
    subprocess.Popen(["updateblastdbdaemon.py", "start", "200"])
    today = time.strftime("%Y_%m_%d", time.localtime())
    log_file = os.path.join("/", "home", "primerdesign", "speciesprimer_" + today + ".log")
    if os.path.isfile("/tmp/frontail.pid"):
        with open("/tmp/frontail.pid") as f:
            for line in f:
                pid = line.strip()
        pidint = int(pid)
        if isinstance(pidint, int):
            try:
                os.kill(pidint, signal.SIGTERM)
            except ProcessLookupError:
                pass
        os.remove("/tmp/frontail.pid")
    frontail_cmd = [
            "frontail-linux", "-d", "-n", "20", "--pid-path", "/tmp/frontail.pid", log_file]
    subprocess.Popen(frontail_cmd)

def stop_updatedb():
    subprocess.Popen(["updateblastdbdaemon.py", "stop", "200"])

@app.route('/updatedb', methods=['GET', 'POST'])
def updatedb():
    form = DBForm()
    if form.validate_on_submit():
        if form.submit.data:
            start_updatedb()
            time.sleep(1)
        elif form.stop.data:
            stop_updatedb()
    return render_template('dbdownload.html', title='Control BLAST db download', form=form)

def start_updatedb_ref():
    subprocess.Popen(["updateblastdbdaemon_ref.py", "start", "200"])
    today = time.strftime("%Y_%m_%d", time.localtime())
    log_file = os.path.join("/", "home", "primerdesign", "speciesprimer_" + today + ".log")
    if os.path.isfile("/tmp/frontail.pid"):
        with open("/tmp/frontail.pid") as f:
            for line in f:
                pid = line.strip()
        pidint = int(pid)
        if isinstance(pidint, int):
            try:
                os.kill(pidint, signal.SIGTERM)
            except ProcessLookupError:
                pass
        os.remove("/tmp/frontail.pid")
    frontail_cmd = [
            "frontail-linux", "-d", "-n", "20", "--pid-path", "/tmp/frontail.pid", log_file]
    subprocess.Popen(frontail_cmd)

def stop_updatedb_ref():
    subprocess.Popen(["updateblastdbdaemon_ref.py", "stop", "200"])

@app.route('/updatedb', methods=['GET', 'POST'])
def updatedb_ref():
    form = DBForm()
    if form.validate_on_submit():
        if form.submit.data:
            start_updatedb_ref()
            time.sleep(1)
        elif form.stop.data:
            stop_updatedb_ref()
    return render_template('dbdownload.html', title='Control BLAST db download', form=form)

@app.route('/controlrun', methods=['GET', 'POST'])
def controlrun():
    form = ControlRunForm()
    if form.validate_on_submit():
        if form.submit.data:
            start_pipeline()
            time.sleep(1)
        elif form.stop.data:
            stop_pipeline()
    return render_template('controlrun.html', title='Control runs', form=form)

@app.route('/streamlog', methods=['GET', 'POST'])
def streamlog():
    redirect_url ="http://localhost:9001/"
    return redirect(redirect_url)

def allowed_file(filename):
    ALLOWED_EXTENSIONS = ["json"]
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def get_bacteria_taxids():
    if os.path.isfile(tmp_db_path):
        with open(tmp_db_path) as f:
            for line in f:
                tmp_db = json.loads(line)
    from Bio import Entrez
    Entrez.tool = "SpeciesPrimer pipeline"
    Entrez.email = tmp_db['email']
    searchtaxid = Entrez.esearch(db="taxonomy", term="txid2[orgn]", retmax="500000")
    taxidresult = Entrez.read(searchtaxid)
    taxids = taxidresult["IdList"]
    taxids.sort(key=lambda x: int(x))
    filename = os.path.join(pipe_dir, "dictionaries", "2.txids")
    with open(filename, "w") as f:
        for txid in taxids:
            f.write(txid + "\n")

@app.route('/blastdb', methods=['GET', 'POST'])
def blastdb():
    form = DownloadDB()
    if form.validate_on_submit():
        if form.update_blastdb.data is True:
            return redirect(url_for('updatedb'))
        elif form.update_refprok.data is True:
            return redirect(url_for('updatedb_ref'))
        elif form.update_txids.data is True:
            get_bacteria_taxids()
        elif form.get_blastdb.data is True:
            delete = form.delete.data
            db = form.whichdb.data
            update_dict = {'BLAST_DB':{'delete': delete, 'db': db}}
            update_tmp_db(update_dict)
            return redirect(url_for('dbdownload'))

    return render_template('blastdb.html', title='Download nt database', form=form)

@app.route('/pipelineconfig', methods=['GET', 'POST'])
def pipelineconfig():
    form = PipeConfig()
    if form.validate_on_submit():
        if  form.up_list.data:
            f = form.up_list.data
            filename = secure_filename(f.filename)
            if filename:
                f.save(os.path.join(pipe_dir, "dictionaries", "species_list.txt"))
                flash("saved species list")
            else:
                flash("not a supported file")
            return redirect(url_for('pipelineconfig'))
        elif form.up_abbrev.data:
            f = form.up_abbrev.data
            filename = secure_filename(f.filename)
            if filename:
                f.save(os.path.join(pipe_dir, "dictionaries", "genus_abbrev.csv"))
                flash("saved genus abbreviations")
            else:
                flash("not a supported file")
            return redirect(url_for('pipelineconfig'))
        elif form.up_p3sett.data:
            f = form.up_p3sett.data
            filename = secure_filename(f.filename)
            if filename:
                f.save(os.path.join(pipe_dir, "p3parameters"))
                flash("saved primer3 parameters")
            else:
                flash("not a supported file")
        elif form.up_noblast.data:
            f = form.up_noblast.data
            filename = secure_filename(f.filename)
            if filename:
                f.save(os.path.join(pipe_dir, "NO_Blast", "NO_BLAST.gi"))
                flash("saved NO_BLAST.gi")
            else:
                flash("not a supported file")
        elif form.down_list.data is True:
            return send_file(os.path.join(pipe_dir, "dictionaries", "species_list.txt"), as_attachment=True, attachment_filename="species_list.txt")
        elif form.down_abbrev.data is True:
            return send_file(os.path.join(pipe_dir, "dictionaries", "genus_abbrev.csv"), as_attachment=True, attachment_filename="genus_abbrev.csv")
        elif form.down_p3sett.data is True:
            return send_file(os.path.join(pipe_dir, "p3parameters"), as_attachment=True, attachment_filename="p3parameters")
        elif form.down_noblast.data is True:
            return send_file(os.path.join(pipe_dir, "NO_Blast", "NO_BLAST.gi"), as_attachment=True, attachment_filename='NO_BLAST.gi')

        elif form.reset_list.data is True:
            old_file = os.path.join(pipe_dir, "dictionaries", "species_list.txt")
            default = os.path.join(pipe_dir, "default", "species_list.txt")
            shutil.copy(default, old_file)
            flash("Reset species_list.txt")
        elif form.reset_abbrev.data is True:
            old_file = os.path.join(pipe_dir, "dictionaries", "genus_abbrev.csv")
            default = os.path.join(pipe_dir, "default", "genus_abbrev.csv")
            shutil.copy(default, old_file)
            flash("Reset genus_abbrev.csv")
        elif form.reset_noblast.data is True:
            old_file = os.path.join(pipe_dir, "NO_Blast", "NO_BLAST.gi")
            default = os.path.join(pipe_dir, "default", "NO_BLAST.gi")
            shutil.copy(default, old_file)
            flash("Reset NO_BLAST.gi")
        elif form.reset_p3sett.data is True:
            old_file = os.path.join(pipe_dir, "p3parameters")
            default = os.path.join(pipe_dir, "default", "p3parameters")
            shutil.copy(default, old_file)
            flash("Reset p3parameters")

    return render_template('pipelineconfig.html', title='SpeciesPrimer configuration', form=form)

@app.route('/help')
def helppage():
     return render_template('help.html', title='SpeciesPrimer help')
@app.route('/help_pipelineconfig')
def helppipesetup():
    return render_template('help_pipelinesetup.html', title='SpeciesPrimer help')
@app.route('/help_primerdesign')
def helpprimerdesign():
    return render_template('help_primerdesign.html', title='SpeciesPrimer help')
@app.route('/troubleshooting')
def helptroubleshooting():
    return render_template('help_troubleshooting.html', title='SpeciesPrimer help')
@app.route('/docker')
def helpdocker():
    return render_template('help_dockertroubleshooting.html', title='SpeciesPrimer help')
@app.route('/dockersetup')
def helpdockerproxy():
    return render_template('help_dockerproxy.html', title='SpeciesPrimer help')
