import subprocess
from tkinter import CENTER

##### Schrodinger script

#sch_runtime = '/apps/schrodinger2020-3/run' # schrodinger server
confx_runtime = '/apps/schrodinger2020-3/confgenx' # schrodinger server
#sch_runtime = '/apps/schrodinger2021-2/run'
#confx_runtime = '/apps/schrodinger2021-2/confgenx'
path_sch_script = './sch/'
sch_runtime = '/opt/schrodinger/suites2022-2/run'

def procsch_vsquery(input):
    script_path = path_sch_script + 'phypo.py'
    process = subprocess.Popen([sch_runtime, script_path, input], stdout=subprocess.PIPE)
    stdout = process.communicate()[0]
    res = stdout.decode('utf-8')
    # print(res)
    return res


def procsch_featuregen(input):
    script_path = path_sch_script + 'map_features.py'
    process = subprocess.Popen([sch_runtime, script_path, input], stdout=subprocess.PIPE)
    stdout = process.communicate()[0]
    res = stdout.decode('utf-8')
    return res


def procsch_confgenx(input):
    process = subprocess.Popen([confx_runtime, '-m 256', '-S 0', '-WAIT', '-SAVE', '-NOJOBID', input],
                               stdout=subprocess.PIPE).wait(timeout=None)


def procsch_smitomol(input, output):
    script_path = path_sch_script + 'smitomol.py'
    process = subprocess.Popen([sch_runtime, script_path, input, output], stdout=subprocess.PIPE).wait(
        timeout=None)


def procsch_maegztosdf(input, output):
    script_path = path_sch_script + 'maetosdf.py'
    process = subprocess.Popen([sch_runtime, script_path, input, output],
                               stdout=subprocess.PIPE).wait(
        timeout=None)
                               

# pdb -> add Hydrogen -> receptor.mae
def procsch_addh_rec(input, output):
    script_path = path_sch_script + 'protein_addH.py'
    process = subprocess.Popen([sch_runtime, script_path, input, output],
                               stdout=subprocess.PIPE).wait(
        timeout=None)
# receptor.mae -> get -site_center
def procsch_recsite(input, asl):#, output):
    script_path = path_sch_script + 'get_sitecenter.py'
    process = subprocess.Popen([sch_runtime, script_path, input, asl],
                               stdout=subprocess.PIPE)
    stdout = process.communicate()[0]
    site = stdout.decode('utf-8')
    return site

# make recptor pharmacophore
def procsch_recpharm(rec, lig, pdb):#, output):
    cmd = "/opt/schrodinger/suites2022-2/utilities/phase_complex"+" "+rec+" "+lig+" "+"Hypothesis_"+pdb+" "+ "-HR -HRdist 10.0 -RH -RHdist 10.0 -ADdist 10 -HHdist 10 -NPdist 10 -RRdist 10 -xvol  -rhypo"
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    process.wait()
    stdout = process.communicate()
    print(stdout)
    
# protein preparation
def protein_prep(input, output):#, output):
    cmd = "/opt/schrodinger/suites2022-2/utilities/prepwizard" +" "+input+" "+output+" "+"-fillsidechains -disulfides -rehtreat -max_states 1 -epik_pH 7.4 -epik_pHt 2.0 -antibody_cdr_scheme Kabat -samplewater -propka_pH 7.4 -fix -f S-OPLS -rmsd 0.3 -keepfarwat"
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True).wait(timeout=None)

# make glide input asl
def procsch_glide_asl(input, asl):#, output):
    script_path = path_sch_script + 'glide_asl.py'
    process = subprocess.Popen([sch_runtime, script_path, input, asl],
                               stdout=subprocess.PIPE)
    stdout = process.communicate()[0]
    lig_asl = stdout.decode('utf-8')
    return lig_asl
# make grid
def make_grid(input):#, output):
    cmd = "/opt/schrodinger/suites2022-2/glide "+input
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True).wait(timeout=None)

# start glide docking
def glide_start(input):#, output):
    cmd = "/opt/schrodinger/suites2022-2/glide /Users/song-inhyeok/CODING/PROTOTYPE/BIChem/"+input+" -OVERWRITE -adjust -HOST localhost:4"
    #cmd = "/opt/schrodinger/suites2022-2/glide /Users/song-inhyeok/CODING/PROTOTYPE/BIChem/"+input+" -OVERWRITE -adjust -HOST localhost:4"
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True).wait(timeout=None)
#glide pv -> sdf
def procsch_pvsdf(input, output):
    script_path = path_sch_script + 'glide_resultsplit.py'
    process = subprocess.Popen([sch_runtime, script_path, input, output],
                               stdout=subprocess.PIPE).wait(
        timeout=None)
                               
#MD_make_model 
def make_model(input):
    cmd ="/opt/schrodinger/suites2022-1/utilities/multisim -JOBNAME "+ input+ "_setup -m "+ input+ ".msj " +"/Users/song-inhyeok/CODING/PROTOTYPE/BIChem/"+input+"_prep.maegz -o "+input+"-out.cms -maxjob 1 -HOST localhost"
    #cmd ="$SCHRODINGER/utilities/multisim -JOBNAME " + input +"-m" +input + ".msj" +input+ ".mae "+"-o"+ input+".cms -HOST localhost"
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True).wait(timeout=None)

#MD_start
def start_md(input):
    cmd ='/opt/schrodinger/suites2022-1/utilities/multisim  -JOBNAME ' +input+'_MD -HOST gnode2 -maxjob 1 -cpu 1 -m /Users/song-inhyeok/CODING/PROTOTYPE/BIChem/'+input+'_md.msj -c /Users/song-inhyeok/CODING/PROTOTYPE/BIChem/MD_templete.cfg -description "Molecular Dynamics" /Users/song-inhyeok/CODING/PROTOTYPE/BIChem/'+input+'-out.cms -mode umbrella -o '+input+'_MD-out.cms -lic DESMOND_GPGPU:16'
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True).wait(timeout=None)
