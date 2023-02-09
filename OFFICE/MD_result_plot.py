def PL_RMSD(file):
    import pandas as pd
    df = pd.read_fwf(file)
    df['frame']= df['frame#']
    df= df.drop(['frame#',"#"], axis=1)
    p_max = df['Prot_CA'].max()
    p_min = df['Prot_CA'].min()
    l_min = df['Lig_wrt_Protein'].min()
    l_max = df['Lig_wrt_Protein'].max()
    file_name = file.split(".")[0]
    
    import plotly.graph_objects as go
        
    # import make_subplots function from plotly.subplots
    # to make grid of plots
    from plotly.subplots import make_subplots

    # use specs parameter in make_subplots function
    # to create secondary y-axis
    fig = make_subplots(specs=[[{"secondary_y": True}]])

    # plot a scatter chart by specifying the x and y values
    # Use add_trace function to specify secondary_y axes.
    fig.add_trace(
        go.Scatter(x=df["frame"], y=df['Prot_CA'], name="Protein"),
        secondary_y=False)

    # Use add_trace function and specify secondary_y axes = True.
    fig.add_trace(
        go.Scatter(x=df["frame"], y=df['Lig_wrt_Protein'], name="Ligand"),
        secondary_y=True,)

    # Adding title text to the figure
    fig.update_layout(
        title_text="<b>Protein-Ligand RMSD</b>", title_x=0.5
    )

    # Naming x-axis
    fig.update_xaxes(title_text="Frame", range=[-10,1010])

    # Naming y-axes
    fig.update_yaxes(title_text="<b>Protein RMSD (Å)</b>", range=[p_min, p_max], secondary_y=False)
    fig.update_yaxes(title_text="<b>Lignad RMSD (Å)</b>",range=[l_min, l_max], secondary_y=True)

    fig.update_layout(
        autosize=False,
        width=800,
        height=600,)
    #fig.show()
    fig.write_image(f"{file_name}.svg")
    
    
def Pro_RMSF(file):
    import pandas as pd
    df =pd.read_fwf(file)
    df=df.drop('#', axis=1)
    file_name = file.split(".")[0]
    
    max = df['CA'].max()
    min = df['CA'].min()
    import plotly.graph_objects as go
        
    # import make_subplots function from plotly.subplots
    # to make grid of plots
    from plotly.subplots import make_subplots

    # use specs parameter in make_subplots function
    # to create secondary y-axis
    fig = make_subplots(specs=[[{"secondary_y": True}]])

    # plot a scatter chart by specifying the x and y values
    # Use add_trace function to specify secondary_y axes.
    fig.add_trace(
        go.Scatter(x=df["Residue#"], y=df['CA']),
        secondary_y=False)

    # Use add_trace function and specify secondary_y axes = True.

    # Adding title text to the figure
    fig.update_layout(
        title_text="<b>Protein RMSF</b>", title_x=0.5
    )

    # Naming x-axis
    fig.update_xaxes(title_text="Residue Index")

    # Naming y-axes
    fig.update_yaxes(title_text="<b>RMSF (Å)</b>", range=[min, max], secondary_y=False)
    

    fig.update_layout(
        autosize=False,
        width=800,
        height=600,)
    #fig.show()
    fig.write_image(f"{file_name}.svg")
    
def PL_contact(file1,file2, file3, file4):
    import pandas as pd
    
    df1 =pd.read_fwf(file1)
    df1['Frame'] = df1['Frame#']
    df1=df1.drop(["#","Frame#","AtomName","LigandFragment", "LigandAtomName"], axis=1)
    cat_HB = ['H-Bonds' for j in range(1) for i in range(len(df1.index))]
    df1['cat']=cat_HB
    
    df2 =pd.read_fwf(file2)
    df2['Frame'] = df2['Frame#']
    df2=df2.drop(["#","Frame#","LigandFragment"], axis=1)
    cat_HY = ['Hydrophobic' for j in range(1) for i in range(len(df2.index))]
    df2['cat']=cat_HY

    df3 =pd.read_fwf(file3)
    df3['Frame'] = df3['Frame#']
    df3=df3.drop(["#","Frame#","AtomName","LigandFragment", "LigandAtom","Distance"], axis=1)
    cat_ION = ['Ionic' for j in range(1) for i in range(len(df3.index))]
    df3['cat']=cat_ION
    
    df4 =pd.read_fwf(file4)
    df4['Frame'] = df4['Frame#']
    df4=df4.drop(["#","Frame#","AtomName","LigandFragment", "LigandAtom"], axis=1)
    cat_WATER = ['Water bridges' for j in range(1) for i in range(len(df4.index))]
    df4['cat']=cat_WATER
    
    total_df = pd.concat([df1, df2, df3, df4]) 
    
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
    
    
    from plotly.subplots import make_subplots
    import plotly.graph_objects as go
    import plotly.express as px
    from itertools import cycle
    
    fig = px.histogram(final_df, 
            x=final_df['Res'], 
            y=final_df['percent'], 
            color=final_df['Interaction'], 
            barmode="stack"
            )
    fig.update_layout(title="<b>Protein-Ligand Contacts</b>", title_x=0.5)
    fig.update_layout(yaxis_title="Interaction Fraction")
    fig.update_xaxes(type='category', title_text="Resudue_index")
    fig.update_layout(
            autosize=False,
            width=800,
            height=600,)
    fig.write_image("PL-Contact.svg")  
    
    
    
def ENEplot(file):
    import pandas as pd
    df = pd.read_fwf(file,skiprows=9)
    df["9:T (K)"] = df['(K)']
    df['Frame'] = df['0:time (ps)']*0.0125
    df = df.drop(['9:T', '(K)', '#', '0:time (ps)'], axis=1)
    file_name = file.split(".")[0]
    
    
    max = df['1:E   (kcal/mol)'].max()
    min = df['2:E_p (kcal/mol)'].min()

    import plotly.graph_objects as go
    
    # import make_subplots function from plotly.subplots
    # to make grid of plots
    from plotly.subplots import make_subplots
    
    # use specs parameter in make_subplots function
    # to create secondary y-axis
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    
    # plot a scatter chart by specifying the x and y values
    # Use add_trace function to specify secondary_y axes.
    fig.add_trace(
        go.Scatter(x=df["Frame"], y=df['1:E   (kcal/mol)'], name="Total energy"),
        secondary_y=False)
    
    # Use add_trace function and specify secondary_y axes = True.
    fig.add_trace(
        go.Scatter(x=df["Frame"], y=df['2:E_p (kcal/mol)'], name="Poetential energy"),
        secondary_y=True,)
    
    # Adding title text to the figure
    fig.update_layout(
        title_text="<b>Energy_profile</b>", title_x=0.5
    )
    
    # Naming x-axis
    fig.update_xaxes(title_text="Frame")
    
    # Naming y-axes
    fig.update_yaxes(title_text="<b>Total energy</b>", range=[min-2000, max+2000], secondary_y=False)
    fig.update_yaxes(title_text="<b>Potential energy</b>",range=[min-2000, max+2000], secondary_y=True)

    fig.update_layout(
        autosize=False,
        width=800,
        height=600,)
    #fig.show()
    fig.write_image(f"{file_name}.svg")  
