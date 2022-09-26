# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 11:12:49 2022

@author: adm
"""


# Introduction


# Needed packages

from PyBioMed import Pyprotein
import pandas as pd
from pandas import DataFrame
from Bio import SeqIO
import numpy as np
import streamlit as st
from pathlib import Path
import base64
from PIL import Image
import io
import plotly.graph_objects as go

# Streamlit config

#%%

#---------------------------------#
# Page layout
## Page expands to full width
st.set_page_config(page_title='LIDeB Tools - Druggability', page_icon="📏", layout='wide')

######
# Function to put a picture as header   
def img_to_bytes(img_path):
    img_bytes = Path(img_path).read_bytes()
    encoded = base64.b64encode(img_bytes).decode()
    return encoded

image = Image.open('cropped-header.png')
st.image(image)
#####
st.write("[![Website](https://img.shields.io/badge/website-LIDeB-blue)](https://lideb.biol.unlp.edu.ar)[![Twitter Follow](https://img.shields.io/twitter/follow/LIDeB_UNLP?style=social)](https://twitter.com/intent/follow?screen_name=LIDeB_UNLP)")
st.markdown("###### Made in 🐍 and [![this is an image link](https://i.imgur.com/iIOA6kU.png)](https://www.streamlit.io/) with :heart: by [LIDeB](https://twitter.com/LIDeB_UNLP)")
#####
st.subheader(":pushpin:" "About Us")
st.markdown("We are a drug discovery team with an interest in the development of publicly available open-source customizable cheminformatics tools to be used in computer-assisted drug discovery. We belong to the Laboratory of Bioactive Research and Development (LIDeB) of the National University of La Plata (UNLP), Argentina. Our research group is focused on computer-guided drug repurposing and rational discovery of new drug candidates to treat epilepsy and neglected tropical diseases.")


# Introduction
#---------------------------------#

st.write("""
# LIDeB Tools - FaDrA

**It is a free web-application for Druggability Prediction**

Druggability refers to the ability of a given protein to bind with high affinity to small, drug-like molecules.
Assessing the druggability of potential pharmacological targets is of crucial importance before starting a 
target-focused drug discovery campaign. Fast Druggability Assessment (FaDrA) is a druggability prediction web application based 
on four linear classifiers, capable of discriminating druggable from non-druggable targets from complete proteomes in a few minutes,
with acceptable accuracy, based only on the protein sequence.


The tool uses the following packages: [PyBioMed](https://github.com/gadsbyfly/PyBioMed/), [Biopython](https://biopython.org/), [Plotly](https://plotly.com/)

The next workflow summarizes the steps performed by this method:
    
    
""")


image = Image.open('workflow_druggability.png')
st.image(image, caption='Druggability Workflow')


########### OPTIONS #######
# SIDEBAR

# Loading file
st.sidebar.header('Upload your fasta')

archivo = st.sidebar.file_uploader("Upload a fasta file with sequences", type=["fasta"])

applicability_domain_h = st.sidebar.slider('h* for applicability domain', 2, 3, 3, 1)


# Druggability

def druggability_app(archivo):
    seq=[]
    proteinas = []
    for record in SeqIO.parse(archivo,'fasta'):
        seq.append(record.id)
        sequence = str(record.seq)
        sequence = sequence.strip()
        proteinas.append(sequence)
    
    df=pd.DataFrame()
    df2=pd.DataFrame()
    valores=pd.DataFrame()
    a=Pyprotein.PyProtein(str(proteinas[0]))
    a=a.GetCTD()
    b=a.keys()

    for i in b:
    	claves = pd.DataFrame(columns=[i])
    	df = pd.concat([df,claves],sort=False)
    
    for i in proteinas:
    	valor = pd.DataFrame(columns=[i])
    	df2 = pd.concat([df2,valor],sort=False)
    
    df2=df2.T
    
    for i in range(len(proteinas)):
    	protein_class = Pyprotein.PyProtein(str(proteinas[i]))
    	CTD = protein_class.GetCTD()
    	valor = pd.DataFrame(CTD.values())
    	valor = valor.T
    	valores=pd.concat([valores,valor])
    	valores=valores.rename(index={0:str(seq[i])})
    
    b=list(b)
    for i in range(len(b)):
    	valores=valores.rename(columns={i:str(b[i])})
        
    valores=DataFrame(valores[['_SolventAccessibilityD1100','_SecondaryStrD1025','_ChargeD1100','_SolventAccessibilityD1050','_HydrophobicityD3025',
    			  '_SecondaryStrD3001','_ChargeD2075','_PolarityD1001','_PolarizabilityD2050','_ChargeD3100',
                              '_SolventAccessibilityD1100','_ChargeD2075','_NormalizedVDWVD3025','_NormalizedVDWVD3075','_PolarizabilityD3100',
                              '_SolventAccessibilityD1100','_ChargeD2050','_NormalizedVDWVD3075','_SolventAccessibilityD1075','_PolarizabilityD3025']])
    coeficientes=[0.160,0.038,-0.027,-0.037,0.025,
                  -0.109,-0.096,-0.089,0.028,-0.019,
                  0.174,-0.091,0.026,-0.034,0.044,
                  0.160,-0.054,-0.027,-0.040,0.016]
    descriptores=[]
    cuentas=[]
    
    for j in range(len(valores.columns)):
    	v=valores.iloc[:,j].tolist()
    	descriptores.append(v)
    	v=[numero*coeficientes[j] for numero in v]
    	cuentas.append(v)
    
    M1=[]
    M2=[]
    M3=[]
    M4=[]
    
    z=0
    while z<len(valores.index):
    	m1=cuentas[0][z]+cuentas[1][z]+cuentas[2][z]+cuentas[3][z]+cuentas[4][z]-12.232
    	m2=cuentas[5][z]+cuentas[6][z]+cuentas[7][z]+cuentas[8][z]+cuentas[9][z]+8.400
    	m3=cuentas[10][z]+cuentas[11][z]+cuentas[12][z]+cuentas[13][z]+cuentas[14][z]-12.245
    	m4=cuentas[15][z]+cuentas[16][z]+cuentas[17][z]+cuentas[18][z]+cuentas[19][z]-8.019
    	
    	M1.append(m1)
    	M2.append(m2)
    	M3.append(m3)
    	M4.append(m4)
    	z=z+1
    
    corte=0.5
    scores=(M1,M2,M3,M4)
    # st.write(scores)
    predicciones=[] #Lista de listas, cada elemento de la lista es uno de los modelos. Cada uno de los modelos tiene las train predicciones binarizadas.
    for models in scores:
      for j in models:#Transformo a valor de clase a cada una de las predicciones de los modelos.
          float(j)
          if j<corte:
            j='Non-druggable'
            predicciones.append(j)
          else: 
            j='Druggable'
            predicciones.append(j)	
    
    predicciones=[predicciones[j:j+len(valores.index)] for j in range (0,len(predicciones),len(valores.index))]#Separo la lista original en (n°modelos) listas
    
    resultados=DataFrame(valores.index,columns=['Protein_ID'])
    l=1
    for x in predicciones:
    	g=DataFrame(x,columns=['PREDICT_' + str(l)])
    	resultados=pd.concat([resultados,g],axis=1)
    	l=l+1
    
    train=pd.read_csv('train_todo.csv',sep=',')
    
    train=DataFrame(train[['_SolventAccessibilityD1100','_SecondaryStrD1025','_ChargeD1100','_SolventAccessibilityD1050','_HydrophobicityD3025',
    			  '_SecondaryStrD3001','_ChargeD2075','_PolarityD1001','_PolarizabilityD2050','_ChargeD3100',
                              '_SolventAccessibilityD1100','_ChargeD2075','_NormalizedVDWVD3025','_NormalizedVDWVD3075','_PolarizabilityD3100',
                              '_SolventAccessibilityD1100','_ChargeD2050','_NormalizedVDWVD3075','_SolventAccessibilityD1075','_PolarizabilityD3025']])
    
    h=[]
    
    for i in range(4):
    	m=train.iloc[:,i*5:i*5+5]
    	x=m.to_numpy()
    	xt=np.transpose(x)
    	y=np.matmul(xt,x)
    	z=np.linalg.inv(y)
    	for f in range(len(valores.index)):
    		m_=valores.iloc[f,i*5:i*5+5]
    		x_=m_.to_numpy()
    		xt_=np.transpose(x_)
    		y_=np.matmul(x_,z)
    		z_=np.matmul(y_,xt_)
    		hi=z_
    		h.append(hi)
    if applicability_domain_h == 2:
        h2=10/156 # 2k/n
    else:
        h2=15/156 # 3k/n
    
    dominio=[]
    
    for g in h:
        if g<h2:
            w='YES'
            dominio.append(w)
        else:
            w='NO'
            dominio.append(w)

    
    dominio=[dominio[j:j+len(valores.index)] for j in range (0,len(predicciones*len(valores.index)),len(valores.index))]#Separo la lista original en (n°modelos) listas
    v=0
    dominio=DataFrame(dominio).T
    for y in predicciones:
    	dominio=dominio.rename(columns={v:'DA_'+str(v+1)})
    	v=v+1
    resultados_final = pd.concat([resultados,dominio],axis=1)
    resultados_final['name'] = seq
    columnas=["Protein_ID","PREDICT_1","DA_1","PREDICT_2","DA_2","PREDICT_3","DA_3","PREDICT_4","DA_4"]
    resultados_final_crudo = resultados_final[columnas]
    
    scores1 = pd.DataFrame(scores).T
    
    return resultados_final_crudo, scores1


#%%


def druggability_analisis(resultados_final):
    resultados_final1 = resultados_final.copy()
    
    st.markdown(str(resultados_final1.shape[0]) + " proteins were loaded")
    
    resultados_final1.set_index('Protein_ID',inplace=True)
    
    models = [["PREDICT_1","DA_1"],["PREDICT_2","DA_2"],["PREDICT_3","DA_3"],["PREDICT_4","DA_4"]]
    
    data_4_4 = resultados_final1.copy()
    
    for index, row in data_4_4.iterrows():
        druggables_da = []
        druggables_no_da = []
        non_druggables_da = []
        non_druggables_no_da = []

        for model in models:
            pred = row[model[0]]
            da = row[model[1]]
            if pred == "Druggable" and da == "YES":
                druggables_da.append(1)
            if pred == "Druggable" and da == "NO":
                druggables_no_da.append(1)  
            if pred == "Non-druggable" and da == "YES":
                non_druggables_da.append(1)
            if pred == "Non-druggable" and da == "NO":
                 non_druggables_no_da.append(1)               
            else:
                pass
            
        # Si al menos un modelo la predice "NO Drogable" y entra en DA, queda "No druggable"
        if len(non_druggables_da) >= 1:
            data_4_4.loc[index,'CLASSIFICATION'] = "Non-Druggable"  
        
        # Si hay 3 o más modelos que dan fuera del DA, da "Non-conclusive"
        elif (len(druggables_no_da) + len(non_druggables_no_da)) >= 3:
            data_4_4.loc[index,'CLASSIFICATION'] = "Non-Conclusive"  
        
        # Si las 4 dan drogables y está en da, se considera drogable.
        elif len(druggables_da) == 4:
            data_4_4.loc[index,'CLASSIFICATION'] = "Druggable"
        
        # si 2 dan drogables en DA, y dos no drogables fuera del DA, se considera drogable.
        elif len(druggables_da) == 2 and len(non_druggables_no_da) == 2:
            data_4_4.loc[index,'CLASSIFICATION'] = "Druggable"
        else:
            data_4_4.loc[index,'CLASSIFICATION'] = "Druggable"
        
    all_druggable = data_4_4.loc[data_4_4['CLASSIFICATION'] == "Druggable"]
    all_no_druggable = data_4_4.loc[data_4_4['CLASSIFICATION'] == "Non-Druggable"]
    inconclusive_all = data_4_4.loc[data_4_4['CLASSIFICATION'] == "Non-Conclusive"]

    data_pie_total = [all_druggable.shape[0],all_no_druggable.shape[0],inconclusive_all.shape[0]]


    keys =["Druggable", "Non-Druggable","Non-Conclusive"]
    
    # Create plot:
    
    fig = go.Figure(data=[go.Pie(labels=keys, values=data_pie_total)])
    
    fig.update_layout(plot_bgcolor = 'rgb(256,256,256)',
                            title_font = dict(size=25, family='Calibri', color='black'),
                            font =dict(size=20, family='Calibri'),
                            legend_title_font = dict(size=18, family='Calibri', color='black'),
                            legend_font = dict(size=15, family='Calibri', color='black'))

    st.plotly_chart(fig)
    return data_4_4

#%%
def filedownload(df):
    csv = df.to_csv(index=True,header=True)
    b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
    href = f'<a href="data:file/csv;base64,{b64}" download="Druggability_predictions.csv">Download CSV File with the druggability predictions</a>'
    return href

#%%

if archivo is not None:
    run = st.sidebar.button("PREDICT")
    if run == True:
        file_name_st = archivo.name
        archivo_ok = io.TextIOWrapper(archivo)
        resultados_final_crudo, scores1 = druggability_app(archivo_ok)
        all_to_write = druggability_analisis(resultados_final_crudo)        
        
        st.markdown(":point_down: **Here you can download the raw predictions for 3/4 models**", unsafe_allow_html=True)
        st.markdown(filedownload(all_to_write), unsafe_allow_html=True)


# Example file
else:
    st.info('Awaiting for FASTA file to be uploaded.')
    if st.button('Press to use Example Dataset'):
        archivo = open("example.fasta","rb")
        file_name_st = archivo.name
        archivo_ok = io.TextIOWrapper(archivo)
        resultados_final_crudo, scores1 = druggability_app(archivo_ok)
        all_to_write = druggability_analisis(resultados_final_crudo)        
        st.markdown(":point_down: **Here you can download the raw predictions for 3/4 models**", unsafe_allow_html=True)
        st.markdown(filedownload(all_to_write), unsafe_allow_html=True)


st.markdown("""
         **To cite the application, please reference XXXXXXXXX**
         """)




#Footer edit

footer="""<style>
a:link , a:visited{
color: blue;
background-color: transparent;
text-decoration: underline;
}
a:hover,  a:active {
color: red;
background-color: transparent;
text-decoration: underline;
}
.footer {
position: fixed;
left: 0;
bottom: 0;
width: 100%;
background-color: white;
color: black;
text-align: center;
}
</style>
<div class="footer">
<p>Made in  🐍 and <img style='display: ; ' href="https://streamlit.io" src="https://i.imgur.com/iIOA6kU.png" target="_blank"></img> Developed with ❤️ by <a style='display:; text-align: center;' href="https://lideb.biol.unlp.edu.ar/" target="_blank">LIDeB</a></p>
</div>
"""
st.markdown(footer,unsafe_allow_html=True)

