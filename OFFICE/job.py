from PIL import Image
from io import BytesIO
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
import webdriver_manager
from webdriver_manager.chrome import ChromeDriverManager
import glob
import re

def JobMdResult(request, job_id):
    os.system('sudo docker cp nifty_benz:/mdsrv/server/session_index.json ./media/MD_result')
    check_session_info = "media/MD_result/session_index.json"
    check_df = pd.read_json(check_session_info)
    file_list = os.listdir(f"media/MD_result/{job_id}")
    
    target = f"media/MD_result/{job_id}/*.pdb"
    pdb_list = glob.glob(target)

    if str(job_id) not in check_df['name'].tolist():
        pdb_id=pdb_list[0].split("/")[-1].split("_")[0]
        
        #seleniuim MD_result save   
        chrome_options = webdriver.ChromeOptions()

        chrome_options.add_argument('--headless')

        chrome_options.add_argument('--no-sandbox')

        chrome_options.add_argument('--disable-dev-shm-usage')
        chrome_options.add_argument("--window-size=1920x1080")


        driver = webdriver.Chrome(executable_path="chromedriver",chrome_options=chrome_options)

        driver.implicitly_wait(0.5)
        driver.maximize_window()

        driver.get("http://34.64.245.194:9000")

        driver.find_element("xpath",'//*[@id="app"]/div/div/div[1]/div[3]/div/div/div[1]/button[1]').click()

        time.sleep(1)

        driver.find_element("xpath",'//*[@id="app"]/div/div/div[1]/div[3]/div/div/div[2]/div[1]/div[5]/div/button').click()
        time.sleep(1)


        pdb = driver.find_element("xpath",'//*[@id="app"]/div/div/div[1]/div[3]/div/div/div[2]/div[1]/div[5]/div[3]/div/div/input')
        pdb.send_keys(f"/home/inhyeok.song/Ntrophy/media/MD_result/{job_id}/{pdb_id}_MD_result.pdb")

        time.sleep(1)

        xtc = driver.find_element("xpath",'//*[@id="app"]/div/div/div[1]/div[3]/div/div/div[2]/div[1]/div[5]/div[4]/div/div/input')
        xtc.send_keys(f"/home/inhyeok.song/Ntrophy/media/MD_result/{job_id}/{pdb_id}_MD_result.xtc")

        time.sleep(1)

        driver.find_element("xpath",'//*[@id="app"]/div/div/div[1]/div[3]/div/div/div[2]/div[1]/div[5]/div[5]/div/button').click()

        time.sleep(1)

        driver.find_element("xpath",'//*[@id="app"]/div/div/div[1]/div[3]/div/div/div[1]/button[7]').click()
        time.sleep(1)
        print("spaner")
        driver.find_element("xpath",'//*[@id="app"]/div/div/div[1]/div[4]/div/div/div[4]/div[1]/button').click()
        time.sleep(3)
        print("componet")
        driver.find_element("xpath",'//*[@id="app"]/div/div/div[1]/div[4]/div/div/div[4]/div[3]/div[2]/button[1]').click()
        time.sleep(3)
        print('focus')
        driver.find_element("xpath",'//*[@id="app"]/div/div/div[1]/div[4]/div/div/div[4]/div[3]/div[3]/button[2]').click()
        print('remove')
        time.sleep(1)

        
        driver.find_element("xpath",'//*[@id="app"]/div/div/div[1]/div[1]/div/div[6]/button').click()
        driver.find_element("xpath",'//*[@id="app"]/div/div/div[1]/div[6]/div/div/div[2]/div/button').click()
        
        server_url=driver.find_element("xpath",'//*[@id="app"]/div/div/div[1]/div[6]/div/div/div[2]/div[2]/div/input')
        server_url.clear()
        name=driver.find_element("xpath",'//*[@id="app"]/div/div/div[1]/div[6]/div/div/div[2]/div[3]/div/input')
        description=driver.find_element("xpath",'//*[@id="app"]/div/div/div[1]/div[6]/div/div/div[2]/div[4]/div/input')
        source=driver.find_element("xpath",'//*[@id="app"]/div/div/div[1]/div[6]/div/div/div[2]/div[5]/div/input')

        server_url.send_keys("http://34.64.245.194:9001")
        name.send_keys(job_id)
        description.send_keys(job_id)
        source.send_keys(job_id)

        driver.find_element("xpath",'//*[@id="app"]/div/div/div[1]/div[6]/div/div/div[2]/div[6]/button').click()
        time.sleep(10)
        print('upload')
        driver.close()
        
        #streming server 
        print("start make url")
        os.system('sudo docker cp nifty_benz:/mdsrv/server/session_index.json ./media/MD_result')
        
        session_info = "media/MD_result/session_index.json"
        df = pd.read_json(session_info)
        url_id=df[df['name']==str(job_id)]['id'].tolist()[0]
        
        #
        base_url ="http://34.64.245.194:9000/?session-url=http%3A%2F%2F34.64.245.194%3A9001%2Fget%2Fsession%2F"
        iframe_url = base_url +str(url_id)
        print("end make url")
        
        pl_rmsd_dat = f"media/MD_result/{job_id}/PL_RMSD.dat"
        p_rmsf_dat = f"media/MD_result/{job_id}/P_RMSF.dat"
        contact_dat1 = f"media/MD_result/{job_id}/PL-Contacts_HBond.dat"
        contact_dat2 = f"media/MD_result/{job_id}/PL-Contacts_Hydrophobic.dat"
        contact_dat3 = f"media/MD_result/{job_id}/PL-Contacts_Ionic.dat"
        contact_dat4 = f"media/MD_result/{job_id}/PL-Contacts_WaterBridge.dat"
        ene_dat = f"media/MD_result/{job_id}/{pdb_id}_MD.ene"
        mmgb_dat = f"media/MD_result/{job_id}/mmgbsa_{pdb_id}-prime-out.csv"
        
        PL_RMSD(pl_rmsd_dat,job_id)
        Pro_RMSF(p_rmsf_dat,job_id)
        PL_CONTACT(contact_dat1,contact_dat2, contact_dat3, contact_dat4, job_id)
        ENEplot(ene_dat, job_id)
        MMGB(mmgb_dat, job_id)

        PL_RM = f"media/MD_result/{job_id}/PL_RMSD.json"
        with open(PL_RM, 'r') as plrm:
            PL_RM_json = json.load(plrm)
        
        Pro_RM = f"media/MD_result/{job_id}/Pro_RMSF.json"
        with open(Pro_RM, 'r') as prorm:
            Pro_RM_json = json.load(prorm)
            
            
        PL_CT = f"media/MD_result/{job_id}/PL_contact.json"
        with open(PL_CT, 'r') as plct:
            PL_CT_json = json.load(plct)
            
        ENE = f"media/MD_result/{job_id}/ENEplot.json"
        with open(ENE, 'r') as ene:
            ENE_json = json.load(ene)
            
        MMG = f"media/MD_result/{job_id}/mmgvsa.json"
        with open(MMG, 'r') as mmg:
            MMG_json = json.load(mmg)
        
        content = {}
        content['iframe_url'] = iframe_url
        content['PL_RMSD'] = json.dumps(PL_RM_json,ensure_ascii=False)
        content['Pro_RMSF'] = json.dumps(Pro_RM_json,ensure_ascii=False)
        content['PL_CONTACT'] = json.dumps(PL_CT_json,ensure_ascii=False)
        content['ENE'] = json.dumps(ENE_json,ensure_ascii=False )
        content['MMG'] = json.dumps(MMG_json,ensure_ascii=False )
        
        return render(request, 'job/result/MD/list.html', content)
 
    else: 
        print("start make url")
        os.system('sudo docker cp nifty_benz:/mdsrv/server/session_index.json ./media/MD_result')
        
        session_info = "media/MD_result/session_index.json"
        df = pd.read_json(session_info)
       
        url_id=df[df['name']==str(job_id)]['id'].tolist()[0]
        
        #
        base_url ="http://34.64.245.194:9000/?session-url=http%3A%2F%2F34.64.245.194%3A9001%2Fget%2Fsession%2F"
        iframe_url = base_url +str(url_id)
        print("end make url")
        
        pl_rmsd_dat = f"media/MD_result/{job_id}/PL_RMSD.dat"
        p_rmsf_dat = f"media/MD_result/{job_id}/P_RMSF.dat" 
        contact_dat1 = f"media/MD_result/{job_id}/PL-Contacts_HBond.dat"
        contact_dat2 = f"media/MD_result/{job_id}/PL-Contacts_Hydrophobic.dat"
        contact_dat3 = f"media/MD_result/{job_id}/PL-Contacts_Ionic.dat"
        contact_dat4 = f"media/MD_result/{job_id}/PL-Contacts_WaterBridge.dat"
        ene_dat = f"media/MD_result/{job_id}/NPT_templete.ene"
        mmgb_dat = f"media/MD_result/{job_id}/mmgbsa_1v4s-prime-out.csv"
        #PL_RMSD(pl_rmsd_dat,job_id)
        #Pro_RMSF(p_rmsf_dat,job_id)
        #PL_CONTACT(contact_dat1,contact_dat2, contact_dat3, contact_dat4, job_id)
        #ENEplot(ene_dat, job_id)
      
        PL_RM = f"media/MD_result/{job_id}/PL_RMSD.json"
        with open(PL_RM, 'r') as plrm:
            PL_RM_json = json.load(plrm)
        
        Pro_RM = f"media/MD_result/{job_id}/Pro_RMSF.json"
        with open(Pro_RM, 'r') as prorm:
            Pro_RM_json = json.load(prorm)
            
            
        PL_CT = f"media/MD_result/{job_id}/PL_contact.json"
        with open(PL_CT, 'r') as plct:
            PL_CT_json = json.load(plct)
            
        ENE = f"media/MD_result/{job_id}/ENEplot.json"
        with open(ENE, 'r') as ene:
            ENE_json = json.load(ene)
        
        MMG = f"media/MD_result/{job_id}/mmgvsa.json"
        with open(MMG, 'r') as mmg:
            MMG_json = json.load(mmg)
        
        job_conformer = Job_conformer.objects.filter(job_id=str(job_id)).first()
       
        content = {}
        content['iframe_url'] = iframe_url
        content['PL_RMSD'] = json.dumps(PL_RM_json,ensure_ascii=False)
        content['Pro_RMSF'] = json.dumps(Pro_RM_json,ensure_ascii=False)
        content['PL_CONTACT'] = json.dumps(PL_CT_json,ensure_ascii=False)
        content['ENE'] = json.dumps(ENE_json,ensure_ascii=False )
        content['MMG'] = json.dumps(MMG_json,ensure_ascii=False )
        content['select'] = job_conformer.select
        return render(request, 'job/result/MD/list.html', content)
    
    
# TODO: 로직 구현
@login_required
def MD_analysis(request, job_id):
    
    batch_id=request.GET['batch_id']
    
    output_path = f'./media/MD_result/{job_id}'
    
    createDirectory(output_path)
    
    job  = Job.objects.filter(batch=batch_id).filter(type=9).all()
    
    job_id_list = []
    job_name_info = {}
    file_list = []
    index_list = []

    md_result_info = []    
    root_path = "media/MD_result"
    for jo in job:
        job_id_list.append(jo.id)
        job_name_info[jo.id] = jo.job_name
    job_id_list.sort()

    
    for ji in job_id_list:
        jb_name = job_name_info[ji]
        job_conf =Job_conformer.objects.filter(job_id=str(ji)).first()
        
        smiles=job_conf.conformer.molecule.canonical_smiles
        
        
        mol_3d = job_conf.mol_3d
        mol_3d = mol_3d.lstrip()
        score_info = Job_result.objects.filter(job_conformer=str(job_conf.id)).filter(x="dG Average").first()
        
        mol = Chem.MolFromSmiles(smiles)

        AllChem.Compute2DCoords(mol)

        fig = Draw.MolToImage(mol, size=(200,200))
        img = io.BytesIO()
        fig.save(img, "PNG")
        encoded_img = b64encode(img.getvalue())
        decoded_img = encoded_img.decode('utf-8')
        image_data = f"data:imagepng;base64,{decoded_img}"

        md_result_info.append({'image': image_data, 'score': score_info.value, 'mol_3d': mol_3d, 'name':jb_name})
        
        file_path = root_path +"/" + str(ji) +"/PL_contact.json"
        file_list.append(file_path)
        index_list.append(job_id_list.index(ji) + 1)
        
        with open(file_path, 'r') as file:
            data = json.load(file)
            for ran in range(len(data)):
                data[ran]['id'] = job_id_list.index(ji) + 1
        with open(file_path, 'w') as file:
            json.dump(data, file)    
    
    df = pd.concat(map(pd.read_json,file_list), ignore_index=True)
    df_HB = df.copy()
    df_HY = df.copy()
    df_IN = df.copy()
    df_WA = df.copy()


    df_HB=df_HB.drop(['Hydrophobic', 'Ionic', 'Water bridges'], axis=1)
    df_HY=df_HY.drop(['H-Bonds','Ionic', 'Water bridges'], axis=1)
    df_IN=df_IN.drop(['H-Bonds','Hydrophobic', 'Water bridges'], axis=1)
    df_WA = df_WA.drop(['H-Bonds','Hydrophobic', 'Ionic'], axis=1)
    

    df['total'] = df.apply(add, axis=1)
    df=df.drop(['H-Bonds','Hydrophobic', 'Ionic', 'Water bridges'], axis=1)
    
    total_res_list=df['Res'].tolist()
    del_duple= set(total_res_list)
    del_duple_list = list(del_duple)
    final_df = pd.DataFrame(columns=del_duple_list)
    
    final_df_HB = pd.DataFrame(columns=del_duple_list)
    final_df_HB['job_id']= index_list
    final_df_HB=final_df_HB.set_index('job_id')
    
    final_df_HY = pd.DataFrame(columns=del_duple_list)
    final_df_HY['job_id']= index_list
    final_df_HY=final_df_HY.set_index('job_id')
    
    final_df_IN = pd.DataFrame(columns=del_duple_list)
    final_df_IN['job_id']= index_list
    final_df_IN=final_df_IN.set_index('job_id')
    
    final_df_WA = pd.DataFrame(columns=del_duple_list)
    final_df_WA['job_id']= index_list
    final_df_WA=final_df_WA.set_index('job_id')
    
    
    final_df['job_id']= index_list
    final_df=final_df.set_index('job_id')
    
    for row in range(df.shape[0]):
        ref_res=df.loc[row][0]
        ref_id =df.loc[row][1]
        value = df.loc[row][2]
        final_df.at[ref_id, ref_res] = value
    
    final_df=final_df.fillna(0)
    final_df=final_df.reset_index()
    #final_df.to_json(f"./media/MD_result/{job_id}/compare_sample.json", orient = 'records')
    
    
    copy_test_df =  final_df.copy()
    copy_test_df = copy_test_df.set_index('job_id')
    test_col_list = copy_test_df.columns.tolist()
    
    fig = sns.clustermap(copy_test_df, annot = False, cmap = 'RdYlBu_r', figsize=(15,15))#, metric="correlation")#, vmin = -1, vmax = 1)
    dgram = fig.dendrogram_col.dendrogram

    new_col_list = []
    new_row_list = []

    for col_id in fig.dendrogram_col.reordered_ind:
        new_col_list.append(test_col_list[col_id])


    for row_id in fig.dendrogram_row.reordered_ind:
        new_row_list.append(str(row_id+1))
        
    new_row_list = [int (i) for i in new_row_list]

    #x axis = row

    row_string_list = []
    row_index_list =[]

    for row_idx in new_row_list:
        row_string_list.append(str(row_idx))
        row_index_list.append(new_row_list.index(row_idx))

    row_df = pd.DataFrame()
    row_df['sort'] = row_index_list
    row_df['name'] = row_string_list

    row_info = row_df.to_dict('records')
 

    #col
   
    col_list=new_col_list
    col_index_list = []
    for col in col_list:
        col_index_list.append(col_list.index(col))    
        
    col_info_df = pd.DataFrame()
    col_info_df['sort'] = col_index_list
    col_info_df['name'] = col_list
    col_info=col_info_df.to_dict('records')



    from itertools import product
    total_cobi=list(product(new_row_list, new_col_list[::-1]))
    source_list = []
    target_list = []
    value_list = []
    for t in total_cobi:
        source_list.append(new_row_list.index(t[0]))
        target_list.append(new_col_list.index(t[1]))
        value_list.append(copy_test_df.at[t[0], t[1]])
        
        
    link_info_df = pd.DataFrame()
    link_info_df['source'] = source_list
    link_info_df['target'] = target_list
    link_info_df['value'] = value_list

    link_dict = link_info_df.to_dict('records')

    test_dict = {}
    test_dict['row_nodes']= row_info
    test_dict['col_nodes'] = col_info
    test_dict['links'] = link_dict

    with open(f"media/MD_result/{job_id}/md_analysis.json", "w") as outfile:
        json.dump(test_dict, outfile, indent = 4)
    
    anal_data = f"media/MD_result/{job_id}/md_analysis.json"
    with open(anal_data, 'r') as mdal:
        md_al_json = json.load(mdal)
    
    #hbonds

    for row in range(df_HB.shape[0]):    
        ref_res=df_HB.loc[row][0]
        value = df_HB.loc[row][1]
        ref_id =df_HB.loc[row][2]
        final_df_HB.at[ref_id, ref_res] = value

    final_df_HB=final_df_HB.fillna(0)
    test_col_HB = final_df_HB.columns.tolist()
    
    
    fig_HB = sns.clustermap(final_df_HB, annot = False, cmap = 'RdYlBu_r', figsize=(15,15))#, metric="correlation")#, vmin = -1, vmax = 1)
    dgram = fig_HB.dendrogram_col.dendrogram

    HB_col_list = []
    HB_row_list = []

    for col_id in fig_HB.dendrogram_col.reordered_ind:
        HB_col_list.append(test_col_HB[col_id])


    for row_id in fig_HB.dendrogram_row.reordered_ind:
        HB_row_list.append(str(row_id+1))


    HB_row_list = [int (i) for i in HB_row_list]

    #x axis = row

    row_string_list = []
    row_index_list =[]

    for row_idx in HB_row_list:
        row_string_list.append(str(row_idx))
        row_index_list.append(HB_row_list.index(row_idx))

    row_df = pd.DataFrame()
    row_df['sort'] = row_index_list
    row_df['name'] = row_string_list

    row_info = row_df.to_dict('records')

    #col

    col_list=HB_col_list
    col_index_list = []
    for col in col_list:
        col_index_list.append(col_list.index(col))    
        
    col_info_df = pd.DataFrame()
    col_info_df['sort'] = col_index_list
    col_info_df['name'] = col_list
    col_info=col_info_df.to_dict('records')



    from itertools import product
    total_cobi=list(product(HB_row_list, HB_col_list))
    source_list = []
    target_list = []
    value_list = []
    for t in total_cobi:
        
        source_list.append(HB_row_list.index(t[0]))
        target_list.append(HB_col_list.index(t[1]))
        value_list.append( final_df_HB.at[t[0], t[1]])
        
        
    link_info_df = pd.DataFrame()
    link_info_df['source'] = source_list
    link_info_df['target'] = target_list
    link_info_df['value'] = value_list

    link_dict = link_info_df.to_dict('records')

    test_dict = {}
    test_dict['row_nodes']= row_info
    test_dict['col_nodes'] = col_info
    test_dict['links'] = link_dict


    with open(f"media/MD_result/{job_id}/HB_total.json", "w") as outfile:
        json.dump(test_dict, outfile, indent = 4)
        
    HB_data = f"media/MD_result/{job_id}/HB_total.json"
    with open(HB_data, 'r') as HB:
        md_HB_json = json.load(HB)
 
    
    
    #Hydro

    for row in range(df_HY.shape[0]):    
        ref_res=df_HY.loc[row][0]
        value = df_HY.loc[row][1]
        ref_id =df_HY.loc[row][2]
        final_df_HY.at[ref_id, ref_res] = value

    final_df_HY=final_df_HY.fillna(0)
    test_col_HY = final_df_HY.columns.tolist()
    
    
    fig_HY = sns.clustermap(final_df_HY, annot = False, cmap = 'RdYlBu_r', figsize=(15,15))#, metric="correlation")#, vmin = -1, vmax = 1)
    #fig_HY.savefig(f"media/MD_result/{job_id}/HY.png")
    dgram = fig_HY.dendrogram_col.dendrogram
    
    HY_col_list = []
    HY_row_list = []

    for col_id in fig_HY.dendrogram_col.reordered_ind:
        HY_col_list.append(test_col_HY[col_id])


    for row_id in fig_HY.dendrogram_row.reordered_ind:
        HY_row_list.append(str(row_id+1))


    HY_row_list = [int (i) for i in HY_row_list]

    #x axis = row

    row_string_list = []
    row_index_list =[]

    for row_idx in HY_row_list:
        row_string_list.append(str(row_idx))
        row_index_list.append(HB_row_list.index(row_idx))

    row_df = pd.DataFrame()
    row_df['sort'] = row_index_list
    row_df['name'] = row_string_list

    row_info = row_df.to_dict('records')

    #col

    col_list=HY_col_list
    col_index_list = []
    for col in col_list:
        col_index_list.append(col_list.index(col))    
        
    col_info_df = pd.DataFrame()
    col_info_df['sort'] = col_index_list
    col_info_df['name'] = col_list
    col_info=col_info_df.to_dict('records')



    from itertools import product
    total_cobi=list(product(HY_row_list, HY_col_list))
    source_list = []
    target_list = []
    value_list = []
    for t in total_cobi:
        
        source_list.append(HY_row_list.index(t[0]))
        target_list.append(HY_col_list.index(t[1]))
        value_list.append( final_df_HY.at[t[0], t[1]])
        
        
    link_info_df = pd.DataFrame()
    link_info_df['source'] = source_list
    link_info_df['target'] = target_list
    link_info_df['value'] = value_list

    link_dict = link_info_df.to_dict('records')

    test_dict = {}
    test_dict['row_nodes']= row_info
    test_dict['col_nodes'] = col_info
    test_dict['links'] = link_dict


    with open(f"media/MD_result/{job_id}/HY_total.json", "w") as outfile:
        json.dump(test_dict, outfile, indent = 4)
        
    HY_data = f"media/MD_result/{job_id}/HY_total.json"
    with open(HY_data, 'r') as HYP:
        md_HY_json = json.load(HYP)



    #IONIC
    for row in range(df_IN.shape[0]):    
        ref_res=df_IN.loc[row][0]
        value = df_IN.loc[row][1]
        ref_id =df_IN.loc[row][2]
        final_df_IN.at[ref_id, ref_res] = value

    final_df_IN=final_df_IN.fillna(0)
    test_col_IN = final_df_IN.columns.tolist()
    fig_IN = sns.clustermap(final_df_IN, annot = False, cmap = 'RdYlBu_r', figsize=(15,15))#, metric="correlation")#, vmin = -1, vmax = 1)
    dgram = fig_IN.dendrogram_col.dendrogram

    IN_col_list = []
    IN_row_list = []

    for col_id in fig_IN.dendrogram_col.reordered_ind:
        IN_col_list.append(test_col_IN[col_id])


    for row_id in fig_IN.dendrogram_row.reordered_ind:
        IN_row_list.append(str(row_id+1))


    IN_row_list = [int (i) for i in IN_row_list]

    #x axis = row

    row_string_list = []
    row_index_list =[]

    for row_idx in IN_row_list:
        row_string_list.append(str(row_idx))
        row_index_list.append(IN_row_list.index(row_idx))

    row_df = pd.DataFrame()
    row_df['sort'] = row_index_list
    row_df['name'] = row_string_list

    row_info = row_df.to_dict('records')


    #col

    col_list=IN_col_list
    col_index_list = []
    for col in col_list:
        col_index_list.append(col_list.index(col))    
        
    col_info_df = pd.DataFrame()
    col_info_df['sort'] = col_index_list
    col_info_df['name'] = col_list
    col_info=col_info_df.to_dict('records')



    from itertools import product
    total_cobi=list(product(IN_row_list, IN_col_list[::-1]))
    source_list = []
    target_list = []
    value_list = []
    for t in total_cobi:
        source_list.append(IN_row_list.index(t[0]))
        target_list.append(IN_col_list.index(t[1]))
        value_list.append( final_df_IN.at[t[0], t[1]])
        
        
    link_info_df = pd.DataFrame()
    link_info_df['source'] = source_list
    link_info_df['target'] = target_list
    link_info_df['value'] = value_list

    link_dict = link_info_df.to_dict('records')

    test_dict = {}
    test_dict['row_nodes']= row_info
    test_dict['col_nodes'] = col_info
    test_dict['links'] = link_dict


    with open(f"media/MD_result/{job_id}/IN_total.json", "w") as outfile:
            json.dump(test_dict, outfile, indent = 4)
        
    IN_data = f"media/MD_result/{job_id}/IN_total.json"
    with open(IN_data, 'r') as noic:
        md_IN_json = json.load(noic)

    #WATER
    for row in range(df_WA.shape[0]):    
        ref_res=df_WA.loc[row][0]
        value = df_WA.loc[row][1]
        ref_id =df_WA.loc[row][2]
        final_df_WA.at[ref_id, ref_res] = value

    final_df_WA=final_df_WA.fillna(0)
    test_col_WA = final_df_WA.columns.tolist()
    fig_WA = sns.clustermap(final_df_WA, annot = False, cmap = 'RdYlBu_r', figsize=(15,15))#, metric="correlation")#, vmin = -1, vmax = 1)
    dgram = fig_WA.dendrogram_col.dendrogram

    WA_col_list = []
    WA_row_list = []

    for col_id in fig_WA.dendrogram_col.reordered_ind:
        WA_col_list.append(test_col_WA[col_id])


    for row_id in fig_WA.dendrogram_row.reordered_ind:
        WA_row_list.append(str(row_id+1))


    WA_row_list = [int (i) for i in WA_row_list]

    #x axis = row

    row_string_list = []
    row_index_list =[]

    for row_idx in WA_row_list:
        row_string_list.append(str(row_idx))
        row_index_list.append(WA_row_list.index(row_idx))

    row_df = pd.DataFrame()
    row_df['sort'] = row_index_list
    row_df['name'] = row_string_list

    row_info = row_df.to_dict('records')


    #col

    col_list=WA_col_list
    col_index_list = []
    for col in col_list:
        col_index_list.append(col_list.index(col))    
        
    col_info_df = pd.DataFrame()
    col_info_df['sort'] = col_index_list
    col_info_df['name'] = col_list
    col_info=col_info_df.to_dict('records')



    from itertools import product
    total_cobi=list(product(WA_row_list, WA_col_list))
    source_list = []
    target_list = []
    value_list = []
    for t in total_cobi:
        source_list.append(WA_row_list.index(t[0]))
        target_list.append(WA_col_list.index(t[1]))
        value_list.append( final_df_WA.at[t[0], t[1]])
        
        
    link_info_df = pd.DataFrame()
    link_info_df['source'] = source_list
    link_info_df['target'] = target_list
    link_info_df['value'] = value_list

    link_dict = link_info_df.to_dict('records')

    test_dict = {}
    test_dict['row_nodes']= row_info
    test_dict['col_nodes'] = col_info
    test_dict['links'] = link_dict


    with open(f"media/MD_result/{job_id}/WA_total.json", "w") as outfile:
            json.dump(test_dict, outfile, indent = 4)
        
    WA_data = f"media/MD_result/{job_id}/WA_total.json"
    with open(WA_data, 'r') as WA:
        md_WA_json = json.load(WA)

    res_dis = []
    for res in new_col_list:
        res_dis.append('(:' + "A" + ' and (' + str(res.split("_")[1]) + ')) ')
    
    
    name_HB = []
    for i  in  range(len(md_HB_json['col_nodes'])):
        name_HB.append(md_HB_json['col_nodes'][i]['name'])

    
    name_HY = []
    for i  in  range(len(md_HY_json['col_nodes'])):
        name_HY.append(md_HY_json['col_nodes'][i]['name'])

    name_IN = []
    for i  in  range(len(md_IN_json['col_nodes'])):
        name_IN.append(md_IN_json['col_nodes'][i]['name'])
    
    name_WA = []
    for i  in  range(len(md_WA_json['col_nodes'])):
        name_WA.append(md_WA_json['col_nodes'][i]['name'])
    
  
    res_asl = ",".join(str(e) for e in res_dis)
    content = {}   
    content['analysis_info'] = json.dumps(md_al_json,ensure_ascii=False)
    content['total_res'] = json.dumps(sorted(new_col_list, key=lambda s: int(re.search(r'\d+', s).group())))
    content['md_result'] = json.dumps(md_result_info)
    content['res_dis'] = json.dumps(res_asl)
    content['Hbond'] = json.dumps(md_HB_json, ensure_ascii=False)
    content['Hydro'] = json.dumps(md_HY_json, ensure_ascii=False)
    content['Ion'] = json.dumps(md_IN_json, ensure_ascii=False)
    content['WB'] = json.dumps(md_WA_json, ensure_ascii=False)
    content['HB_col'] =json.dumps(sorted(name_HB, key=lambda s: int(re.search(r'\d+', s).group())))
    content['HY_col'] =json.dumps(sorted(name_HY, key=lambda s: int(re.search(r'\d+', s).group())))
    content['IN_col'] =json.dumps(sorted(name_IN, key=lambda s: int(re.search(r'\d+', s).group())))
    content['WA_col'] =json.dumps(sorted(name_WA, key=lambda s: int(re.search(r'\d+', s).group())))
    
    return render(request, 'job/result/MD/analysis.html', content)
    
   
def add(row):
    return row[1]+row[2]+row[3]+row[4]

# TODO: 로직 구현
@login_required
def retro_sys(request):
    mol_id=request.POST['mol_id']
      
    gce_instance = 'root@nt-chemdb-slurm-cluster-controller'
    gce_project = 'nettargets-chemdb-v1'
    gce_zone = 'asia-northeast3-b'
    gce_query_path = '/home/NCAP_test'
    # gce_output_path = '/home/Blast/bastp.out'
    output_path = './media'

    #cmd_scp = "gcloud compute ssh %s:%s --zone %s --project %s" % (gce_instance, gce_query_path, gce_zone, gce_project)
    #process = subprocess.Popen(cmd_scp.split(), stdout=subprocess.PIPE).wait(timeout=None)
    
    gce_proc_cd = 'cd /home/NCAP_test/ &&'
    gce_proc_py = '/home/hyungsik_jo/anaconda3/envs/upload/bin/python'
    gce_proc_file = '/home/NCAP_test/chemdb.py'
    options = f'--cmd retro --pdb {str(mol_id)}'
    
    cmd_screening = "--command=%s %s %s %s" % (gce_proc_cd, gce_proc_py, gce_proc_file, options)
    process = subprocess.Popen(
        ['gcloud', 'compute', 'ssh', gce_instance, '--zone=asia-northeast3-b', cmd_screening], stdout=subprocess.PIPE).wait()
    
    gce_output_path = '/home/NCAP_test/route_file.tar.gz'
    #gce_output_path = '/home/BIChem/ligand/ligand_pharm.txt'

    cmd_out_scp = "gcloud compute scp %s:%s %s --zone %s --project %s" % (gce_instance, gce_output_path, output_path, gce_zone, gce_project)
    process = subprocess.Popen(cmd_out_scp.split(), stdout=subprocess.PIPE).wait(timeout=None)
    os.system(f'mkdir ./media/retro_result/{str(mol_id)}')
    os.system(f'tar -xvf media/route_file.tar.gz -C ./media/retro_result/{str(mol_id)}')
    os.remove('media/route_file.tar.gz')
    file_list = os.listdir(f"./media/retro_result/{str(mol_id)}/root/data/retro_result/") 
    
    data = []
    for file in file_list:
        full_name  = f"./media/retro_result/{str(mol_id)}/root/data/retro_result/"+file
        data.append(full_name)  
    
    img_source =[]    
    for da in data:
        
        out = BytesIO()

        with Image.open(da) as img:
            img.save(out, format="png")
    
            encoded_img = b64encode(out.getvalue())
            decoded_img = encoded_img.decode('utf-8')
            image_data = f"data:imagepng;base64,{decoded_img}"
            img_source.append(image_data)
    
    content ={'img':img_source }
    return JsonResponse(content)

def MD_ResultSave(request):

    prep_info = dict(request.POST.lists())
    start_info=prep_info["info[0][]"]
    end_info=prep_info["info[1][]"]
    ave_info=prep_info["info[2][]"]
    std_info=prep_info["info[3][]"]
    job_id = prep_info["info[4][]"][0]
    
    std_info[0] = 'dG SD'

    result_list = [start_info, end_info, ave_info, std_info]
    ## result_list example = [['start', 10, 'ns'], ['end', 20, 'ns'], ['dG average', -51.3515, 'kcal/mol(MMGBSA)'], ['dG standard deviation', 3.5084, 'kcal/mol(MMGBSA)']]
    
    job_conformer = Job_conformer.objects.filter(job_id=str(job_id)).first()
    for result in result_list:   
        column_name = result[0]    
        value = result[1]
        unit_name = result[2]
        unit = Unit.objects.filter(title = unit_name).first()
        print(column_name, value, unit_name, unit)
        
        job_result = Job_result()
        
        job_result.job_conformer_id = job_conformer.id
        job_result.index = 0
        job_result.x = column_name
        job_result.value = value
        job_result.unit = unit
        
        job_result.save()
    
    job_conformer.select = True
    job_conformer.save()

@csrf_exempt
def MD_Analysis_Redraw(request):
    anal_info = dict(request.POST.lists())

    batch_id=anal_info['batch'][0]
    an_job = anal_info['job[]'] 
    an_job_list = []
    for aj in an_job:
        an_job_list.append(int(aj))
    an_col_list = anal_info['col[]']
    
    
    #output_path = f'./media/MD_result/{job_id}'
    
    #createDirectory(output_path)
    
    job  = Job.objects.filter(batch=batch_id).filter(type=9).all()
    job_id_list = []
    file_list = []
    index_list = []

    md_result_info = []    
    root_path = "media/MD_result"
    for jo in job:
        job_id_list.append(jo.id)
        
    job_id_list.sort()
    
    
    for ji in job_id_list: 
        file_path = root_path +"/" + str(ji) +"/PL_contact.json"
        file_list.append(file_path)
        index_list.append(job_id_list.index(ji) + 1)
    
    df = pd.concat(map(pd.read_json,file_list), ignore_index=True)
    df['total'] = df.apply(add, axis=1)
    df=df.drop(['H-Bonds','Hydrophobic', 'Ionic', 'Water bridges'], axis=1)
    
    total_res_list=df['Res'].tolist()
    del_duple= set(total_res_list)
    del_duple_list = list(del_duple)
    final_df = pd.DataFrame(columns=del_duple_list)
    
    final_df['job_id']= index_list
    final_df=final_df.set_index('job_id')
    
    
    for row in range(df.shape[0]):
        ref_res=df.loc[row][0]
        ref_id =df.loc[row][1]
        value = df.loc[row][2]
        final_df.at[ref_id, ref_res] = value
    
    final_df=final_df.fillna(0)
   
    
    drop_index_list = list(set(index_list) - set(an_job_list))
    drop_col_list = list(set(del_duple_list) - set(an_col_list))
 
    final_df.drop(drop_index_list, inplace=True)
    final_df.drop(drop_col_list, axis=1, inplace=True)
    final_df=final_df.reset_index()
    new_job_list = final_df['job_id'].tolist()
    
    copy_test_df =  final_df.copy()
    copy_test_df = copy_test_df.set_index('job_id')

    test_col_list = copy_test_df.columns.tolist()
    fig = sns.clustermap(copy_test_df, annot = False, cmap = 'RdYlBu_r', figsize=(15,15))#, metric="correlation")#, vmin = -1, vmax = 1)
    dgram = fig.dendrogram_col.dendrogram

    new_col_list = []
    new_row_list = []

    for col_id in fig.dendrogram_col.reordered_ind:
        new_col_list.append(test_col_list[col_id])


    for row_id in fig.dendrogram_row.reordered_ind:
        new_row_list.append(new_job_list[row_id])
        
    
    
    new_row_list = [int (i) for i in new_row_list]

    #x axis = row

    row_string_list = []
    row_index_list =[]
    
    for row_idx in new_row_list:
        row_string_list.append(str(row_idx))
        row_index_list.append(new_row_list.index(row_idx))
    
   
    row_df = pd.DataFrame()
    row_df['sort'] = row_index_list
    row_df['name'] = row_string_list

    row_info = row_df.to_dict('records')
    
    #col
   
    col_list=new_col_list
   
    col_index_list = []
    for col in col_list:
        col_index_list.append(col_list.index(col))    
        
    col_info_df = pd.DataFrame()
    col_info_df['sort'] = col_index_list
    col_info_df['name'] = col_list
    col_info=col_info_df.to_dict('records')

    

    from itertools import product
    total_cobi=list(product(new_row_list, new_col_list[::-1]))
    source_list = []
    target_list = []
    value_list = []


    for t in total_cobi:
        source_list.append(new_row_list.index(t[0]))
        target_list.append(new_col_list.index(t[1]))
        value_list.append(copy_test_df.at[t[0], t[1]])
        
        
    link_info_df = pd.DataFrame()
    link_info_df['source'] = source_list
    link_info_df['target'] = target_list
    link_info_df['value'] = value_list

    link_dict = link_info_df.to_dict('records')

    test_dict = {}
    test_dict['row_nodes']= row_info
    test_dict['col_nodes'] = col_info
    test_dict['links'] = link_dict
  
  
    return JsonResponse(test_dict)


@csrf_exempt
def MD_Analysis_Redraw_HB(request):
    anal_info = dict(request.POST.lists())

    batch_id=anal_info['batch'][0]
    an_job = anal_info['job[]'] 
    an_job_list = []
    for aj in an_job:
        an_job_list.append(int(aj))
    an_col_list = anal_info['col[]']
    
    
    #output_path = f'./media/MD_result/{job_id}'
    
    #createDirectory(output_path)
    
    job  = Job.objects.filter(batch=batch_id).filter(type=9).all()
    job_id_list = []
    file_list = []
    index_list = []

    md_result_info = []    
    root_path = "media/MD_result"
    for jo in job:
        job_id_list.append(jo.id)
        
    job_id_list.sort()
    
    
    for ji in job_id_list: 
        file_path = root_path +"/" + str(ji) +"/PL_contact.json"
        file_list.append(file_path)
        index_list.append(job_id_list.index(ji) + 1)
    
    df = pd.concat(map(pd.read_json,file_list), ignore_index=True)
    
    df=df.drop(['Hydrophobic', 'Ionic', 'Water bridges'], axis=1)
    
    total_res_list=df['Res'].tolist()
    del_duple= set(total_res_list)
    del_duple_list = list(del_duple)
    final_df = pd.DataFrame(columns=del_duple_list)
    
    final_df['job_id']= index_list
    final_df=final_df.set_index('job_id')
    
    
    for row in range(df.shape[0]):
        ref_res=df.loc[row][0]
        ref_id =df.loc[row][2]
        value = df.loc[row][1]
        final_df.at[ref_id, ref_res] = value
    
    final_df=final_df.fillna(0)
   
    
    drop_index_list = list(set(index_list) - set(an_job_list))
    drop_col_list = list(set(del_duple_list) - set(an_col_list))
 
    final_df.drop(drop_index_list, inplace=True)
    final_df.drop(drop_col_list, axis=1, inplace=True)
    final_df=final_df.reset_index()
    new_job_list = final_df['job_id'].tolist()
    
    copy_test_df =  final_df.copy()
    copy_test_df = copy_test_df.set_index('job_id')

    test_col_list = copy_test_df.columns.tolist()
    fig = sns.clustermap(copy_test_df, annot = False, cmap = 'RdYlBu_r', figsize=(15,15))#, metric="correlation")#, vmin = -1, vmax = 1)
    dgram = fig.dendrogram_col.dendrogram

    new_col_list = []
    new_row_list = []

    for col_id in fig.dendrogram_col.reordered_ind:
        new_col_list.append(test_col_list[col_id])


    for row_id in fig.dendrogram_row.reordered_ind:
        new_row_list.append(new_job_list[row_id])
        
    
    
    new_row_list = [int (i) for i in new_row_list]

    #x axis = row

    row_string_list = []
    row_index_list =[]
    
    for row_idx in new_row_list:
        row_string_list.append(str(row_idx))
        row_index_list.append(new_row_list.index(row_idx))
    
   
    row_df = pd.DataFrame()
    row_df['sort'] = row_index_list
    row_df['name'] = row_string_list

    row_info = row_df.to_dict('records')
    
    #col
   
    col_list=new_col_list
   
    col_index_list = []
    for col in col_list:
        col_index_list.append(col_list.index(col))    
        
    col_info_df = pd.DataFrame()
    col_info_df['sort'] = col_index_list
    col_info_df['name'] = col_list
    col_info=col_info_df.to_dict('records')

    

    from itertools import product
    total_cobi=list(product(new_row_list, new_col_list))
    source_list = []
    target_list = []
    value_list = []


    for t in total_cobi:
        source_list.append(new_row_list.index(t[0]))
        target_list.append(new_col_list.index(t[1]))
        value_list.append(copy_test_df.at[t[0], t[1]])
        
        
    link_info_df = pd.DataFrame()
    link_info_df['source'] = source_list
    link_info_df['target'] = target_list
    link_info_df['value'] = value_list

    link_dict = link_info_df.to_dict('records')

    test_dict = {}
    test_dict['row_nodes']= row_info
    test_dict['col_nodes'] = col_info
    test_dict['links'] = link_dict
  
  
    return JsonResponse(test_dict)


@csrf_exempt
def MD_Analysis_Redraw_HY(request):
    anal_info = dict(request.POST.lists())

    batch_id=anal_info['batch'][0]
    an_job = anal_info['job[]'] 
    an_job_list = []
    for aj in an_job:
        an_job_list.append(int(aj))
    an_col_list = anal_info['col[]']
    
    
    #output_path = f'./media/MD_result/{job_id}'
    
    #createDirectory(output_path)
    
    job  = Job.objects.filter(batch=batch_id).filter(type=9).all()
    job_id_list = []
    file_list = []
    index_list = []

    md_result_info = []    
    root_path = "media/MD_result"
    for jo in job:
        job_id_list.append(jo.id)
        
    job_id_list.sort()
    
    
    for ji in job_id_list: 
        file_path = root_path +"/" + str(ji) +"/PL_contact.json"
        file_list.append(file_path)
        index_list.append(job_id_list.index(ji) + 1)
    
    df = pd.concat(map(pd.read_json,file_list), ignore_index=True)
    df=df.drop(['H-Bonds','Ionic', 'Water bridges'], axis=1)
    
    total_res_list=df['Res'].tolist()
    del_duple= set(total_res_list)
    del_duple_list = list(del_duple)
    final_df = pd.DataFrame(columns=del_duple_list)
    
    final_df['job_id']= index_list
    final_df=final_df.set_index('job_id')
    
    
    for row in range(df.shape[0]):
        
        ref_res=df.loc[row][0]
        ref_id =df.loc[row][2]
        value = df.loc[row][1]
        final_df.at[ref_id, ref_res] = value
    
    final_df=final_df.fillna(0)
   
    
    drop_index_list = list(set(index_list) - set(an_job_list))
    drop_col_list = list(set(del_duple_list) - set(an_col_list))
 
    final_df.drop(drop_index_list, inplace=True)
    final_df.drop(drop_col_list, axis=1, inplace=True)
    final_df=final_df.reset_index()
    new_job_list = final_df['job_id'].tolist()
    
    copy_test_df =  final_df.copy()
    copy_test_df = copy_test_df.set_index('job_id')

    test_col_list = copy_test_df.columns.tolist()
    fig = sns.clustermap(copy_test_df, annot = False, cmap = 'RdYlBu_r', figsize=(15,15))#, metric="correlation")#, vmin = -1, vmax = 1)
    dgram = fig.dendrogram_col.dendrogram

    new_col_list = []
    new_row_list = []

    for col_id in fig.dendrogram_col.reordered_ind:
        new_col_list.append(test_col_list[col_id])


    for row_id in fig.dendrogram_row.reordered_ind:
        new_row_list.append(new_job_list[row_id])
        
    
    
    new_row_list = [int (i) for i in new_row_list]

    #x axis = row

    row_string_list = []
    row_index_list =[]
    
    for row_idx in new_row_list:
        row_string_list.append(str(row_idx))
        row_index_list.append(new_row_list.index(row_idx))
    
   
    row_df = pd.DataFrame()
    row_df['sort'] = row_index_list
    row_df['name'] = row_string_list

    row_info = row_df.to_dict('records')
    
    #col
   
    col_list=new_col_list
   
    col_index_list = []
    for col in col_list:
        col_index_list.append(col_list.index(col))    
        
    col_info_df = pd.DataFrame()
    col_info_df['sort'] = col_index_list
    col_info_df['name'] = col_list
    col_info=col_info_df.to_dict('records')

    

    from itertools import product
    total_cobi=list(product(new_row_list, new_col_list))
    source_list = []
    target_list = []
    value_list = []


    for t in total_cobi:
        source_list.append(new_row_list.index(t[0]))
        target_list.append(new_col_list.index(t[1]))
        value_list.append(copy_test_df.at[t[0], t[1]])
        
        
    link_info_df = pd.DataFrame()
    link_info_df['source'] = source_list
    link_info_df['target'] = target_list
    link_info_df['value'] = value_list

    link_dict = link_info_df.to_dict('records')

    test_dict = {}
    test_dict['row_nodes']= row_info
    test_dict['col_nodes'] = col_info
    test_dict['links'] = link_dict
  
  
    return JsonResponse(test_dict)


@csrf_exempt
def MD_Analysis_Redraw_IN(request):
    anal_info = dict(request.POST.lists())

    batch_id=anal_info['batch'][0]
    an_job = anal_info['job[]'] 
    an_job_list = []
    for aj in an_job:
        an_job_list.append(int(aj))
    an_col_list = anal_info['col[]']
    
    
    #output_path = f'./media/MD_result/{job_id}'
    
    #createDirectory(output_path)
    
    job  = Job.objects.filter(batch=batch_id).filter(type=9).all()
    job_id_list = []
    file_list = []
    index_list = []

    md_result_info = []    
    root_path = "media/MD_result"
    for jo in job:
        job_id_list.append(jo.id)
        
    job_id_list.sort()
    
    
    for ji in job_id_list: 
        file_path = root_path +"/" + str(ji) +"/PL_contact.json"
        file_list.append(file_path)
        index_list.append(job_id_list.index(ji) + 1)
    
    df = pd.concat(map(pd.read_json,file_list), ignore_index=True)
    df=df.drop(['H-Bonds','Hydrophobic', 'Water bridges'], axis=1)
    
    total_res_list=df['Res'].tolist()
    del_duple= set(total_res_list)
    del_duple_list = list(del_duple)
    final_df = pd.DataFrame(columns=del_duple_list)
    
    final_df['job_id']= index_list
    final_df=final_df.set_index('job_id')
    
    
    for row in range(df.shape[0]):
        ref_res=df.loc[row][0]
        ref_id =df.loc[row][2]
        value = df.loc[row][1]
        final_df.at[ref_id, ref_res] = value
    
    final_df=final_df.fillna(0)
   
    
    drop_index_list = list(set(index_list) - set(an_job_list))
    drop_col_list = list(set(del_duple_list) - set(an_col_list))
 
    final_df.drop(drop_index_list, inplace=True)
    final_df.drop(drop_col_list, axis=1, inplace=True)
    final_df=final_df.reset_index()
    new_job_list = final_df['job_id'].tolist()
    
    copy_test_df =  final_df.copy()
    copy_test_df = copy_test_df.set_index('job_id')

    test_col_list = copy_test_df.columns.tolist()
    fig = sns.clustermap(copy_test_df, annot = False, cmap = 'RdYlBu_r', figsize=(15,15))#, metric="correlation")#, vmin = -1, vmax = 1)
    dgram = fig.dendrogram_col.dendrogram

    new_col_list = []
    new_row_list = []

    for col_id in fig.dendrogram_col.reordered_ind:
        new_col_list.append(test_col_list[col_id])


    for row_id in fig.dendrogram_row.reordered_ind:
        new_row_list.append(new_job_list[row_id])
        
    
    
    new_row_list = [int (i) for i in new_row_list]

    #x axis = row

    row_string_list = []
    row_index_list =[]
    
    for row_idx in new_row_list:
        row_string_list.append(str(row_idx))
        row_index_list.append(new_row_list.index(row_idx))
    
   
    row_df = pd.DataFrame()
    row_df['sort'] = row_index_list
    row_df['name'] = row_string_list

    row_info = row_df.to_dict('records')
    
    #col
   
    col_list=new_col_list
   
    col_index_list = []
    for col in col_list:
        col_index_list.append(col_list.index(col))    
        
    col_info_df = pd.DataFrame()
    col_info_df['sort'] = col_index_list
    col_info_df['name'] = col_list
    col_info=col_info_df.to_dict('records')

    

    from itertools import product
    total_cobi=list(product(new_row_list, new_col_list))
    source_list = []
    target_list = []
    value_list = []


    for t in total_cobi:
        source_list.append(new_row_list.index(t[0]))
        target_list.append(new_col_list.index(t[1]))
        value_list.append(copy_test_df.at[t[0], t[1]])
        
        
    link_info_df = pd.DataFrame()
    link_info_df['source'] = source_list
    link_info_df['target'] = target_list
    link_info_df['value'] = value_list

    link_dict = link_info_df.to_dict('records')

    test_dict = {}
    test_dict['row_nodes']= row_info
    test_dict['col_nodes'] = col_info
    test_dict['links'] = link_dict
  
  
    return JsonResponse(test_dict)

@csrf_exempt
def MD_Analysis_Redraw_WA(request):
    anal_info = dict(request.POST.lists())

    batch_id=anal_info['batch'][0]
    an_job = anal_info['job[]'] 
    an_job_list = []
    for aj in an_job:
        an_job_list.append(int(aj))
    an_col_list = anal_info['col[]']
    
    
    #output_path = f'./media/MD_result/{job_id}'
    
    #createDirectory(output_path)
    
    job  = Job.objects.filter(batch=batch_id).filter(type=9).all()
    job_id_list = []
    file_list = []
    index_list = []

    md_result_info = []    
    root_path = "media/MD_result"
    for jo in job:
        job_id_list.append(jo.id)
        
    job_id_list.sort()
    
    
    for ji in job_id_list: 
        file_path = root_path +"/" + str(ji) +"/PL_contact.json"
        file_list.append(file_path)
        index_list.append(job_id_list.index(ji) + 1)
    
    df = pd.concat(map(pd.read_json,file_list), ignore_index=True)
    df['total'] = df.apply(add, axis=1)
    df=df.drop(['H-Bonds','Hydrophobic', 'Ionic'], axis=1)
    
    total_res_list=df['Res'].tolist()
    del_duple= set(total_res_list)
    del_duple_list = list(del_duple)
    final_df = pd.DataFrame(columns=del_duple_list)
    
    final_df['job_id']= index_list
    final_df=final_df.set_index('job_id')
    
    
    for row in range(df.shape[0]):
        ref_res=df.loc[row][0]
        ref_id =df.loc[row][2]
        value = df.loc[row][1]
        final_df.at[ref_id, ref_res] = value
    
    final_df=final_df.fillna(0)
   
    
    drop_index_list = list(set(index_list) - set(an_job_list))
    drop_col_list = list(set(del_duple_list) - set(an_col_list))
 
    final_df.drop(drop_index_list, inplace=True)
    final_df.drop(drop_col_list, axis=1, inplace=True)
    final_df=final_df.reset_index()
    new_job_list = final_df['job_id'].tolist()
    
    copy_test_df =  final_df.copy()
    copy_test_df = copy_test_df.set_index('job_id')

    test_col_list = copy_test_df.columns.tolist()
    fig = sns.clustermap(copy_test_df, annot = False, cmap = 'RdYlBu_r', figsize=(15,15))#, metric="correlation")#, vmin = -1, vmax = 1)
    dgram = fig.dendrogram_col.dendrogram

    new_col_list = []
    new_row_list = []

    for col_id in fig.dendrogram_col.reordered_ind:
        new_col_list.append(test_col_list[col_id])


    for row_id in fig.dendrogram_row.reordered_ind:
        new_row_list.append(new_job_list[row_id])
        
    
    
    new_row_list = [int (i) for i in new_row_list]

    #x axis = row

    row_string_list = []
    row_index_list =[]
    
    for row_idx in new_row_list:
        row_string_list.append(str(row_idx))
        row_index_list.append(new_row_list.index(row_idx))
    
   
    row_df = pd.DataFrame()
    row_df['sort'] = row_index_list
    row_df['name'] = row_string_list

    row_info = row_df.to_dict('records')
    
    #col
   
    col_list=new_col_list
   
    col_index_list = []
    for col in col_list:
        col_index_list.append(col_list.index(col))    
        
    col_info_df = pd.DataFrame()
    col_info_df['sort'] = col_index_list
    col_info_df['name'] = col_list
    col_info=col_info_df.to_dict('records')

    

    from itertools import product
    total_cobi=list(product(new_row_list, new_col_list[::-1]))
    source_list = []
    target_list = []
    value_list = []


    for t in total_cobi:
        source_list.append(new_row_list.index(t[0]))
        target_list.append(new_col_list.index(t[1]))
        value_list.append(copy_test_df.at[t[0], t[1]])
        
        
    link_info_df = pd.DataFrame()
    link_info_df['source'] = source_list
    link_info_df['target'] = target_list
    link_info_df['value'] = value_list

    link_dict = link_info_df.to_dict('records')

    test_dict = {}
    test_dict['row_nodes']= row_info
    test_dict['col_nodes'] = col_info
    test_dict['links'] = link_dict
  
  
    return JsonResponse(test_dict)


def PL_RMSD(file, job_id):
    import pandas as pd
    df = pd.read_fwf(file)
    df['Time']= df['frame#']*0.08
    df= df.drop(['frame#',"#", 'Prot_Backbone', 'Prot_Sidechain','Prot_All_Heavy', 'Lig_wrt_Ligand'], axis=1)
    #df.to_csv(f"./media/MD_result/{job_id}/PL_RMSD.csv", index=False)
    df.to_json(f"./media/MD_result/{job_id}/PL_RMSD.json", orient = 'records')
    
def Pro_RMSF(file, job_id):
    df =pd.read_fwf(file)
    df['Residue'] = df['Residue#']
    df=df.drop(['#', 'Residue#','Chain', 'ResName',  'LigandContact', 'Backbone', 'Sidechain','All_Heavy','B-factor'], axis=1)
    #df.to_csv(f"./media/MD_result/{job_id}/Pro_RMSF.csv", index=False)
    df.to_json(f"./media/MD_result/{job_id}/Pro_RMSF.json", orient = 'records')
    
def PL_CONTACT(file,file1,file2,file3, job_id):
    df1 =pd.read_fwf(file)
    df1['Frame'] = df1['Frame#']
    df1=df1.drop(["#","Frame#","AtomName","LigandFragment", "LigandAtomName"], axis=1)
    cat_HB = ['H-Bonds' for j in range(1) for i in range(len(df1.index))]
    df1['cat']=cat_HB

    df2 =pd.read_fwf(file1)
    df2['Frame'] = df2['Frame#']
    df2=df2.drop(["#","Frame#","LigandFragment"], axis=1)
    cat_HY = ['Hydrophobic' for j in range(1) for i in range(len(df2.index))]
    df2['cat']=cat_HY

    df3 =pd.read_fwf(file2)
    df3['Frame'] = df3['Frame#']
    df3=df3.drop(["#","Frame#","AtomName","LigandFragment", "LigandAtom","Distance"], axis=1)
    cat_ION = ['Ionic' for j in range(1) for i in range(len(df3.index))]
    df3['cat']=cat_ION

    df4 =pd.read_fwf(file3)
    df4['Frame'] = df4['Frame#']
    df4=df4.drop(["#","Frame#","AtomName","LigandFragment", "LigandAtom"], axis=1)
    cat_WATER = ['Water bridges' for j in range(1) for i in range(len(df4.index))]
    df4['cat']=cat_WATER

        
    concat_list = []
    #1111
    if df1.shape[0]!=0 and df2.shape[0]!=0 and df3.shape[0]!=0 and df2.shape[0]!=0:
        concat_list.append(df1)
        concat_list.append(df2) 
        concat_list.append(df3) 
        concat_list.append(df4)  

        total_df = pd.concat(concat_list) 

        name = []

        for h,z in zip(total_df['ResName'],total_df['Residue#']):
            res_name = h+"_"+str(z)
            name.append(res_name)

        total_df['Res']=name
        total_df= total_df.drop(["Residue#","ResName"], axis=1)

        group_df =pd.DataFrame(total_df.groupby("Res")['cat'].value_counts())

        group_df['count']=group_df['cat']
        group_df=group_df.drop('cat', axis=1)
        final_df=group_df.T


        import numpy as np

        change_list=final_df.columns.tolist()

        per_list =[] 
        for change in change_list:   
            value=final_df[change[0]][change[1]].values
            if change[1]=="H-Bonds":
                per = value /(total_df['cat'].value_counts()[0])
                fi_per=np.round(per,4)
            elif change[1]=="Water bridges":
                per = value /(total_df['cat'].value_counts()[1])
                fi_per=np.round(per,4)
            elif change[1]=="Hydrophobic":
                per = value /(total_df['cat'].value_counts()[2])
                fi_per=np.round(per,4)
            else :
                per = value /(total_df['cat'].value_counts()[3])
                fi_per=np.round(per,4)
            per_list.append(fi_per[0])
            
        final_df.loc['percent']=per_list

        final_df=final_df.drop('count')
        final_df=final_df.T

        final_df=final_df.reset_index()
        final_df['Interaction'] =final_df['cat'] 
        final_df= final_df.drop('cat', axis=1)
        copy_df=final_df.pivot(index='Res', columns='Interaction', values='percent')
        copy_df = copy_df.fillna(0)
        copy_df=copy_df.reset_index()
        copy_df.to_json(f"./media/MD_result/{job_id}/PL_contact.json", orient = 'records')

    #1011
    elif df1.shape[0]!=0 and df2.shape[0]==0 and df3.shape[0]!=0 and df2.shape[0]!=0:
        concat_list.append(df1)
        concat_list.append(df3) 
        concat_list.append(df4)  
        total_df = pd.concat(concat_list) 

        name = []

        for h,z in zip(total_df['ResName'],total_df['Residue#']):
            res_name = h+"_"+str(z)
            name.append(res_name)

        total_df['Res']=name
        total_df= total_df.drop(["Residue#","ResName"], axis=1)

        group_df =pd.DataFrame(total_df.groupby("Res")['cat'].value_counts())

        group_df['count']=group_df['cat']
        group_df=group_df.drop('cat', axis=1)
        final_df=group_df.T


        import numpy as np

        change_list=final_df.columns.tolist()

        per_list =[] 
        for change in change_list:   
            value=final_df[change[0]][change[1]].values
            if change[1]=="H-Bonds":
                per = value /(total_df['cat'].value_counts()[0])
                fi_per=np.round(per,4)
            elif change[1]=="Water bridges":
                per = value /(total_df['cat'].value_counts()[1])
                fi_per=np.round(per,4)
            elif change[1]=="Hydrophobic":
                per = 0
                fi_per=np.round(per,4)
            else :
                per = value /(total_df['cat'].value_counts()[2])
            per_list.append(fi_per[0])
            
        final_df.loc['percent']=per_list

        final_df=final_df.drop('count')
        final_df=final_df.T

        final_df=final_df.reset_index()
        final_df['Interaction'] =final_df['cat'] 
        final_df= final_df.drop('cat', axis=1)
        copy_df=final_df.pivot(index='Res', columns='Interaction', values='percent')
        copy_df = copy_df.fillna(0)
        copy_df=copy_df.reset_index()
        copy_df.to_json(f"./media/MD_result/{job_id}/PL_contact.json", orient = 'records')
    #1101    
    elif df1.shape[0]!=0 and df2.shape[0]!=0 and df3.shape[0]==0 and df2.shape[0]!=0:
        concat_list.append(df1)
        concat_list.append(df2) 
        concat_list.append(df4)     
        total_df = pd.concat(concat_list) 

        name = []

        for h,z in zip(total_df['ResName'],total_df['Residue#']):
            res_name = h+"_"+str(z)
            name.append(res_name)

        total_df['Res']=name
        total_df= total_df.drop(["Residue#","ResName"], axis=1)

        group_df =pd.DataFrame(total_df.groupby("Res")['cat'].value_counts())

        group_df['count']=group_df['cat']
        group_df=group_df.drop('cat', axis=1)
        final_df=group_df.T


        import numpy as np

        change_list=final_df.columns.tolist()

        per_list =[] 
        for change in change_list:   
            value=final_df[change[0]][change[1]].values
            if change[1]=="H-Bonds":
                per = value /(total_df['cat'].value_counts()[0])
                fi_per=np.round(per,4)
            elif change[1]=="Water bridges":
                per = value /(total_df['cat'].value_counts()[1])
                fi_per=np.round(per,4)
            elif change[1]=="Hydrophobic":
                per = value /(total_df['cat'].value_counts()[2])
                fi_per=np.round(per,4)
            else :
                per = 0
                fi_per=np.round(per,4)
            per_list.append(fi_per[0])
            
        final_df.loc['percent']=per_list

        final_df=final_df.drop('count')
        final_df=final_df.T

        final_df=final_df.reset_index()
        final_df['Interaction'] =final_df['cat'] 
        final_df= final_df.drop('cat', axis=1)
        copy_df=final_df.pivot(index='Res', columns='Interaction', values='percent')
        copy_df = copy_df.fillna(0)
        copy_df=copy_df.reset_index()
        copy_df.to_json(f"./media/MD_result/{job_id}/PL_contact.json", orient = 'records')

    #1110
    elif df1.shape[0]!=0 and df2.shape[0]!=0 and df3.shape[0]!=0 and df2.shape[0]==0:   
        concat_list.append(df1)
        concat_list.append(df2) 
        concat_list.append(df3) 
        total_df = pd.concat(concat_list) 

        name = []

        for h,z in zip(total_df['ResName'],total_df['Residue#']):
            res_name = h+"_"+str(z)
            name.append(res_name)

        total_df['Res']=name
        total_df= total_df.drop(["Residue#","ResName"], axis=1)

        group_df =pd.DataFrame(total_df.groupby("Res")['cat'].value_counts())

        group_df['count']=group_df['cat']
        group_df=group_df.drop('cat', axis=1)
        final_df=group_df.T


        import numpy as np

        change_list=final_df.columns.tolist()

        per_list =[] 
        for change in change_list:   
            value=final_df[change[0]][change[1]].values
            if change[1]=="H-Bonds":
                per = value /(total_df['cat'].value_counts()[0])
                fi_per=np.round(per,4)
            elif change[1]=="Water bridges":
                per = 0
                fi_per=np.round(per,4)
            elif change[1]=="Hydrophobic":
                per = value /(total_df['cat'].value_counts()[1])
                fi_per=np.round(per,4)
            else :
                per = value /(total_df['cat'].value_counts()[2])
                fi_per=np.round(per,4)
            per_list.append(fi_per[0])
            
        final_df.loc['percent']=per_list

        final_df=final_df.drop('count')
        final_df=final_df.T

        final_df=final_df.reset_index()
        final_df['Interaction'] =final_df['cat'] 
        final_df= final_df.drop('cat', axis=1)
        copy_df=final_df.pivot(index='Res', columns='Interaction', values='percent')
        copy_df = copy_df.fillna(0)
        copy_df=copy_df.reset_index()
        copy_df.to_json(f"./media/MD_result/{job_id}/PL_contact.json", orient = 'records')
    #0111
    elif df1.shape[0]==0 and df2.shape[0]!=0 and df3.shape[0]!=0 and df2.shape[0]!=0:   
        concat_list.append(df2)
        concat_list.append(df3) 
        concat_list.append(df4) 
        total_df = pd.concat(concat_list) 

        name = []

        for h,z in zip(total_df['ResName'],total_df['Residue#']):
            res_name = h+"_"+str(z)
            name.append(res_name)

        total_df['Res']=name
        total_df= total_df.drop(["Residue#","ResName"], axis=1)

        group_df =pd.DataFrame(total_df.groupby("Res")['cat'].value_counts())

        group_df['count']=group_df['cat']
        group_df=group_df.drop('cat', axis=1)
        final_df=group_df.T


        import numpy as np

        change_list=final_df.columns.tolist()

        per_list =[] 
        for change in change_list:   
            value=final_df[change[0]][change[1]].values
            if change[1]=="H-Bonds":
                per = 0
                fi_per=np.round(per,4)
            elif change[1]=="Water bridges":
                per = value /(total_df['cat'].value_counts()[0])
                fi_per=np.round(per,4)
            elif change[1]=="Hydrophobic":
                per = value /(total_df['cat'].value_counts()[1])
                fi_per=np.round(per,4)
            else :
                per = value /(total_df['cat'].value_counts()[2])
                fi_per=np.round(per,4)
            per_list.append(fi_per[0])
            
        final_df.loc['percent']=per_list

        final_df=final_df.drop('count')
        final_df=final_df.T

        final_df=final_df.reset_index()
        final_df['Interaction'] =final_df['cat'] 
        final_df= final_df.drop('cat', axis=1)
        copy_df=final_df.pivot(index='Res', columns='Interaction', values='percent')
        copy_df = copy_df.fillna(0)
        copy_df=copy_df.reset_index()
        copy_df.to_json(f"./media/MD_result/{job_id}/PL_contact.json", orient = 'records')
    #1100
    elif df1.shape[0]!=0 and df2.shape[0]!=0 and df3.shape[0]==0 and df2.shape[0]==0:   
        concat_list.append(df1)
        concat_list.append(df2) 
        total_df = pd.concat(concat_list) 

        name = []

        for h,z in zip(total_df['ResName'],total_df['Residue#']):
            res_name = h+"_"+str(z)
            name.append(res_name)

        total_df['Res']=name
        total_df= total_df.drop(["Residue#","ResName"], axis=1)

        group_df =pd.DataFrame(total_df.groupby("Res")['cat'].value_counts())

        group_df['count']=group_df['cat']
        group_df=group_df.drop('cat', axis=1)
        final_df=group_df.T


        import numpy as np

        change_list=final_df.columns.tolist()

        per_list =[] 
        for change in change_list:   
            value=final_df[change[0]][change[1]].values
            if change[1]=="H-Bonds":
                per = value /(total_df['cat'].value_counts()[0])
                fi_per=np.round(per,4)
            elif change[1]=="Water bridges":
                per = value /(total_df['cat'].value_counts()[2])
                fi_per=np.round(per,4)
            elif change[1]=="Hydrophobic":
                per = 0
                fi_per=np.round(per,4)
            else :
                per = 0
                fi_per=np.round(per,4)
            per_list.append(fi_per[0])
            
        final_df.loc['percent']=per_list

        final_df=final_df.drop('count')
        final_df=final_df.T

        final_df=final_df.reset_index()
        final_df['Interaction'] =final_df['cat'] 
        final_df= final_df.drop('cat', axis=1)
        copy_df=final_df.pivot(index='Res', columns='Interaction', values='percent')
        copy_df = copy_df.fillna(0)
        copy_df=copy_df.reset_index()
        copy_df.to_json(f"./media/MD_result/{job_id}/PL_contact.json", orient = 'records')
    #1010
    elif df1.shape[0]!=0 and df2.shape[0]==0 and df3.shape[0]!=0 and df2.shape[0]==0:   
        concat_list.append(df1)
        concat_list.append(df3) 
        total_df = pd.concat(concat_list) 

        name = []

        for h,z in zip(total_df['ResName'],total_df['Residue#']):
            res_name = h+"_"+str(z)
            name.append(res_name)

        total_df['Res']=name
        total_df= total_df.drop(["Residue#","ResName"], axis=1)

        group_df =pd.DataFrame(total_df.groupby("Res")['cat'].value_counts())

        group_df['count']=group_df['cat']
        group_df=group_df.drop('cat', axis=1)
        final_df=group_df.T


        import numpy as np

        change_list=final_df.columns.tolist()

        per_list =[] 
        for change in change_list:   
            value=final_df[change[0]][change[1]].values
            if change[1]=="H-Bonds":
                per = value /(total_df['cat'].value_counts()[0])
                fi_per=np.round(per,4)
            elif change[1]=="Water bridges":
                per = 0
                fi_per=np.round(per,4)
            elif change[1]=="Hydrophobic":
                per = 0
                fi_per=np.round(per,4)
            else :
                per = value /(total_df['cat'].value_counts()[2])
                fi_per=np.round(per,4)
            per_list.append(fi_per[0])
            
        final_df.loc['percent']=per_list

        final_df=final_df.drop('count')
        final_df=final_df.T

        final_df=final_df.reset_index()
        final_df['Interaction'] =final_df['cat'] 
        final_df= final_df.drop('cat', axis=1)
        copy_df=final_df.pivot(index='Res', columns='Interaction', values='percent')
        copy_df = copy_df.fillna(0)
        copy_df=copy_df.reset_index()
        copy_df.to_json(f"./media/MD_result/{job_id}/PL_contact.json", orient = 'records')

    #1001
    elif df1.shape[0]!=0 and df2.shape[0]==0 and df3.shape[0]==0 and df2.shape[0]!=0:   
        concat_list.append(df1)
        concat_list.append(df4) 
        total_df = pd.concat(concat_list) 

        name = []

        for h,z in zip(total_df['ResName'],total_df['Residue#']):
            res_name = h+"_"+str(z)
            name.append(res_name)

        total_df['Res']=name
        total_df= total_df.drop(["Residue#","ResName"], axis=1)

        group_df =pd.DataFrame(total_df.groupby("Res")['cat'].value_counts())

        group_df['count']=group_df['cat']
        group_df=group_df.drop('cat', axis=1)
        final_df=group_df.T


        import numpy as np

        change_list=final_df.columns.tolist()

        per_list =[] 
        for change in change_list:   
            value=final_df[change[0]][change[1]].values
            if change[1]=="H-Bonds":
                per = value /(total_df['cat'].value_counts()[0])
                fi_per=np.round(per,4)
            elif change[1]=="Water bridges":
                per = value /(total_df['cat'].value_counts()[1])
                fi_per=np.round(per,4)
            elif change[1]=="Hydrophobic":
                per = 0
                fi_per=np.round(per,4)
            else :
                per = 0
                fi_per=np.round(per,4)
            per_list.append(fi_per[0])
            
        final_df.loc['percent']=per_list

        final_df=final_df.drop('count')
        final_df=final_df.T

        final_df=final_df.reset_index()
        final_df['Interaction'] =final_df['cat'] 
        final_df= final_df.drop('cat', axis=1)
        copy_df=final_df.pivot(index='Res', columns='Interaction', values='percent')
        copy_df = copy_df.fillna(0)
        copy_df=copy_df.reset_index()
        copy_df.to_json(f"./media/MD_result/{job_id}/PL_contact.json", orient = 'records')

    #0110
    elif df1.shape[0]==0 and df2.shape[0]!=0 and df3.shape[0]!=0 and df2.shape[0]==0:
        concat_list.append(df2)
        concat_list.append(df3)  

        total_df = pd.concat(concat_list) 

        name = []

        for h,z in zip(total_df['ResName'],total_df['Residue#']):
            res_name = h+"_"+str(z)
            name.append(res_name)

        total_df['Res']=name
        total_df= total_df.drop(["Residue#","ResName"], axis=1)
    
        group_df =pd.DataFrame(total_df.groupby("Res")['cat'].value_counts())

        group_df['count']=group_df['cat']
        group_df=group_df.drop('cat', axis=1)
        final_df=group_df.T


        import numpy as np

        change_list=final_df.columns.tolist()

        per_list =[] 
        for change in change_list:   
            value=final_df[change[0]][change[1]].values
            
            if change[1]=="H-Bonds":
                per = 0
                fi_per=np.round(per,4)
            elif change[1]=="Water bridges":
                per = 0
                fi_per=np.round(per,4)
            elif change[1]=="Hydrophobic":
                per = value /(total_df['cat'].value_counts()[0])
                fi_per=np.round(per,4)
            else :
                per = value /(total_df['cat'].value_counts()[1])
                fi_per=np.round(per,4)
            per_list.append(fi_per[0])
        final_df.loc['percent']=per_list

        final_df=final_df.drop('count')
        final_df=final_df.T

        final_df=final_df.reset_index()
        final_df['Interaction'] =final_df['cat'] 
        final_df= final_df.drop('cat', axis=1)
        copy_df=final_df.pivot(index='Res', columns='Interaction', values='percent')
        copy_df = copy_df.fillna(0)
        copy_df=copy_df.reset_index()

        #copy_df.to_csv(f"./media/MD_result/{job_id}/PL_contact.csv", index=False)
        copy_df.to_json(f"./media/MD_result/{job_id}/PL_contact.json", orient = 'records')
        
    #0101
    elif df1.shape[0]==0 and df2.shape[0]!=0 and df3.shape[0]==0 and df2.shape[0]!=0:
        concat_list.append(df2)
        concat_list.append(df4)  

        total_df = pd.concat(concat_list) 

        name = []

        for h,z in zip(total_df['ResName'],total_df['Residue#']):
            res_name = h+"_"+str(z)
            name.append(res_name)

        total_df['Res']=name
        total_df= total_df.drop(["Residue#","ResName"], axis=1)
    
        group_df =pd.DataFrame(total_df.groupby("Res")['cat'].value_counts())

        group_df['count']=group_df['cat']
        group_df=group_df.drop('cat', axis=1)
        final_df=group_df.T


        import numpy as np

        change_list=final_df.columns.tolist()

        per_list =[] 
        for change in change_list:   
            value=final_df[change[0]][change[1]].values
            
            if change[1]=="H-Bonds":
                per = 0
                fi_per=np.round(per,4)
            elif change[1]=="Water bridges":
                per = value /(total_df['cat'].value_counts()[0])
                fi_per=np.round(per,4)
            elif change[1]=="Hydrophobic":
                per = value /(total_df['cat'].value_counts()[1])
                fi_per=np.round(per,4)
            else :
                per = 0
                fi_per=np.round(per,4)
            per_list.append(fi_per[0])
        final_df.loc['percent']=per_list

        final_df=final_df.drop('count')
        final_df=final_df.T

        final_df=final_df.reset_index()
        final_df['Interaction'] =final_df['cat'] 
        final_df= final_df.drop('cat', axis=1)
        copy_df=final_df.pivot(index='Res', columns='Interaction', values='percent')
        copy_df = copy_df.fillna(0)
        copy_df=copy_df.reset_index()
        print(copy_df)
        #copy_df.to_csv(f"./media/MD_result/{job_id}/PL_contact.csv", index=False)
        copy_df.to_json(f"./media/MD_result/{job_id}/PL_contact.json", orient = 'records')
        
    #0011
    elif df1.shape[0]==0 and df2.shape[0]==0 and df3.shape[0]!=0 and df2.shape[0]!=0:
        concat_list.append(df3)
        concat_list.append(df4)  

        total_df = pd.concat(concat_list) 

        name = []

        for h,z in zip(total_df['ResName'],total_df['Residue#']):
            res_name = h+"_"+str(z)
            name.append(res_name)

        total_df['Res']=name
        total_df= total_df.drop(["Residue#","ResName"], axis=1)
    
        group_df =pd.DataFrame(total_df.groupby("Res")['cat'].value_counts())

        group_df['count']=group_df['cat']
        group_df=group_df.drop('cat', axis=1)
        final_df=group_df.T


        import numpy as np

        change_list=final_df.columns.tolist()

        per_list =[] 
        for change in change_list:   
            value=final_df[change[0]][change[1]].values
            
            if change[1]=="H-Bonds":
                per = 0
                fi_per=np.round(per,4)
            elif change[1]=="Water bridges":
                per = value /(total_df['cat'].value_counts()[0])
                fi_per=np.round(per,4)
            elif change[1]=="Hydrophobic":
                per = 0
                fi_per=np.round(per,4)
            else :
                per = value /(total_df['cat'].value_counts()[1])
                fi_per=np.round(per,4)
            per_list.append(fi_per[0])
        final_df.loc['percent']=per_list

        final_df=final_df.drop('count')
        final_df=final_df.T

        final_df=final_df.reset_index()
        final_df['Interaction'] =final_df['cat'] 
        final_df= final_df.drop('cat', axis=1)
        copy_df=final_df.pivot(index='Res', columns='Interaction', values='percent')
        copy_df = copy_df.fillna(0)
        copy_df=copy_df.reset_index()
        print(copy_df)
        #copy_df.to_csv(f"./media/MD_result/{job_id}/PL_contact.csv", index=False)
        copy_df.to_json(f"./media/MD_result/{job_id}/PL_contact.json", orient = 'records')
        
    #1000
    elif df1.shape[0]!=0 and df2.shape[0]==0 and df3.shape[0]==0 and df2.shape[0]==0:
        concat_list.append(df1)
        

        total_df = df1

        name = []

        for h,z in zip(total_df['ResName'],total_df['Residue#']):
            res_name = h+"_"+str(z)
            name.append(res_name)

        total_df['Res']=name
        total_df= total_df.drop(["Residue#","ResName"], axis=1)
    
        group_df =pd.DataFrame(total_df.groupby("Res")['cat'].value_counts())

        group_df['count']=group_df['cat']
        group_df=group_df.drop('cat', axis=1)
        final_df=group_df.T


        import numpy as np

        change_list=final_df.columns.tolist()

        per_list =[] 
        for change in change_list:   
            value=final_df[change[0]][change[1]].values
            
            if change[1]=="H-Bonds":
                per = value /(total_df['cat'].value_counts()[0])
                fi_per=np.round(per,4)
            elif change[1]=="Water bridges":
                per = 0
                fi_per=np.round(per,4)
            elif change[1]=="Hydrophobic":
                per = 0
                fi_per=np.round(per,4)
            else :
                per = 0
                fi_per=np.round(per,4)
            per_list.append(fi_per[0])
        final_df.loc['percent']=per_list

        final_df=final_df.drop('count')
        final_df=final_df.T

        final_df=final_df.reset_index()
        final_df['Interaction'] =final_df['cat'] 
        final_df= final_df.drop('cat', axis=1)
        copy_df=final_df.pivot(index='Res', columns='Interaction', values='percent')
        copy_df = copy_df.fillna(0)
        copy_df=copy_df.reset_index()
        print(copy_df)
        #copy_df.to_csv(f"./media/MD_result/{job_id}/PL_contact.csv", index=False)
        copy_df.to_json(f"./media/MD_result/{job_id}/PL_contact.json", orient = 'records')
        
    #0100
    elif df1.shape[0]==0 and df2.shape[0]!=0 and df3.shape[0]==0 and df2.shape[0]==0:
        concat_list.append(df2)
        

        total_df = df2

        name = []

        for h,z in zip(total_df['ResName'],total_df['Residue#']):
            res_name = h+"_"+str(z)
            name.append(res_name)

        total_df['Res']=name
        total_df= total_df.drop(["Residue#","ResName"], axis=1)
    
        group_df =pd.DataFrame(total_df.groupby("Res")['cat'].value_counts())

        group_df['count']=group_df['cat']
        group_df=group_df.drop('cat', axis=1)
        final_df=group_df.T


        import numpy as np

        change_list=final_df.columns.tolist()

        per_list =[] 
        for change in change_list:   
            value=final_df[change[0]][change[1]].values
            
            if change[1]=="H-Bonds":
                per = 0
                fi_per=np.round(per,4)
            elif change[1]=="Water bridges":
                per = 0
                fi_per=np.round(per,4)
            elif change[1]=="Hydrophobic":
                per = value /(total_df['cat'].value_counts()[0])
                fi_per=np.round(per,4)
            else :
                per = 0
                fi_per=np.round(per,4)
            per_list.append(fi_per[0])
        final_df.loc['percent']=per_list

        final_df=final_df.drop('count')
        final_df=final_df.T

        final_df=final_df.reset_index()
        final_df['Interaction'] =final_df['cat'] 
        final_df= final_df.drop('cat', axis=1)
        copy_df=final_df.pivot(index='Res', columns='Interaction', values='percent')
        copy_df = copy_df.fillna(0)
        copy_df=copy_df.reset_index()
        print(copy_df)
        #copy_df.to_csv(f"./media/MD_result/{job_id}/PL_contact.csv", index=False)
        copy_df.to_json(f"./media/MD_result/{job_id}/PL_contact.json", orient = 'records')
        
    #0010
    elif df1.shape[0]==0 and df2.shape[0]==0 and df3.shape[0]!=0 and df2.shape[0]==0:
        concat_list.append(df3)
        

        total_df = df3

        name = []

        for h,z in zip(total_df['ResName'],total_df['Residue#']):
            res_name = h+"_"+str(z)
            name.append(res_name)

        total_df['Res']=name
        total_df= total_df.drop(["Residue#","ResName"], axis=1)
    
        group_df =pd.DataFrame(total_df.groupby("Res")['cat'].value_counts())

        group_df['count']=group_df['cat']
        group_df=group_df.drop('cat', axis=1)
        final_df=group_df.T


        import numpy as np

        change_list=final_df.columns.tolist()

        per_list =[] 
        for change in change_list:   
            value=final_df[change[0]][change[1]].values
            
            if change[1]=="H-Bonds":
                per = 0
                fi_per=np.round(per,4)
            elif change[1]=="Water bridges":
                per = 0
                fi_per=np.round(per,4)
            elif change[1]=="Hydrophobic":
                per = 0
                fi_per=np.round(per,4)
            else :
                per = value /(total_df['cat'].value_counts()[0])
                fi_per=np.round(per,4)
            per_list.append(fi_per[0])
        final_df.loc['percent']=per_list

        final_df=final_df.drop('count')
        final_df=final_df.T

        final_df=final_df.reset_index()
        final_df['Interaction'] =final_df['cat'] 
        final_df= final_df.drop('cat', axis=1)
        copy_df=final_df.pivot(index='Res', columns='Interaction', values='percent')
        copy_df = copy_df.fillna(0)
        copy_df=copy_df.reset_index()
        print(copy_df)
        #copy_df.to_csv(f"./media/MD_result/{job_id}/PL_contact.csv", index=False)
        copy_df.to_json(f"./media/MD_result/{job_id}/PL_contact.json", orient = 'records')
        
    #0001
    elif df1.shape[0]==0 and df2.shape[0]==0 and df3.shape[0]==0 and df2.shape[0]!=0:
        concat_list.append(df4)
        

        total_df = df4

        name = []

        for h,z in zip(total_df['ResName'],total_df['Residue#']):
            res_name = h+"_"+str(z)
            name.append(res_name)

        total_df['Res']=name
        total_df= total_df.drop(["Residue#","ResName"], axis=1)
    
        group_df =pd.DataFrame(total_df.groupby("Res")['cat'].value_counts())

        group_df['count']=group_df['cat']
        group_df=group_df.drop('cat', axis=1)
        final_df=group_df.T


        import numpy as np

        change_list=final_df.columns.tolist()

        per_list =[] 
        for change in change_list:   
            value=final_df[change[0]][change[1]].values
            
            if change[1]=="H-Bonds":
                per = 0
                fi_per=np.round(per,4)
            elif change[1]=="Water bridges":
                per = value /(total_df['cat'].value_counts()[0])
                fi_per=np.round(per,4)
            elif change[1]=="Hydrophobic":
                per = 0
                fi_per=np.round(per,4)
            else :
                per = 0
                fi_per=np.round(per,4)
            per_list.append(fi_per[0])
        final_df.loc['percent']=per_list

        final_df=final_df.drop('count')
        final_df=final_df.T

        final_df=final_df.reset_index()
        final_df['Interaction'] =final_df['cat'] 
        final_df= final_df.drop('cat', axis=1)
        copy_df=final_df.pivot(index='Res', columns='Interaction', values='percent')
        copy_df = copy_df.fillna(0)
        copy_df=copy_df.reset_index()
        print(copy_df)
        #copy_df.to_csv(f"./media/MD_result/{job_id}/PL_contact.csv", index=False)
        copy_df.to_json(f"./media/MD_result/{job_id}/PL_contact.json", orient = 'records')
        
def ENEplot(file, job_id):
    df = pd.read_fwf(file,skiprows=9)
    filter_list = [ (i * 66) for i in range(1000)]
    #filter_list = [ (i * 66) for i in range(250)]
    df_sample = df.loc[filter_list,:]

    df_sample["9:T (K)"] = df_sample['(K)']
    df_sample['Time'] = df_sample['0:time (ps)']*0.001
    df_sample['E'] = df_sample['1:E   (kcal/mol)']
    df_sample['E_p']=df_sample['2:E_p (kcal/mol)']

    df_sample['T'] = df_sample['(K)']
    df_sample['P'] = df_sample['7:P   (bar)']
    df_sample['V'] = df_sample['8:V   (A^3)']
    df_sample = df_sample.drop(['1:E   (kcal/mol)','2:E_p (kcal/mol)','9:T', '(K)', '#', '0:time (ps)','3:E_k (kcal/mol)',
        '4:E_c (kcal/mol)', '5:E_x (kcal/mol)', '6:E_f (kcal/mol)',
        '7:P   (bar)', '8:V   (A^3)', '9:T (K)'], axis=1)

    #df.to_csv(f"./media/MD_result/{job_id}/ENEplot.csv", index=False)
    df_sample.to_json(f"./media/MD_result/{job_id}/ENEplot.json", orient = 'records')

def MMGB(file, job_id):
    df = pd.read_csv(file)
    df['r_psp_MMGBSA_dG_Bind'].to_json(f"./media/MD_result/{job_id}/mmgvsa.json", orient = 'records')

def createDirectory(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print("Error: Failed to create the directory.")
