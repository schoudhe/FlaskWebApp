import pandas
import gseapy
import os, shutil
import glob
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import seaborn as sns

class Model:
    def __init__(self,data,AC_data, biospan_data):
        print("Setting up the DataBase...")

        #self.metadata = pandas.read_csv(data,low_memory=False) didnt work memorywise so uploaded by chunks
        mylist = []
        for chunk in  pandas.read_csv(data, chunksize=20000):
            mylist.append(chunk)
        self.metadata= pandas.concat(mylist, axis= 0)
        del mylist

        self.current_data=self.metadata
        self.AC_data = pandas.read_csv(AC_data)
        self.current_AC = self.AC_data
        self.BS_data  = pandas.read_csv(biospan_data)
        self.current_BS = self.BS_data

    def save_data(self):
        self.current_data.to_csv("output_data.csv", index=False)

    def reset_filter(self):
        self.current_data = self.metadata
        self.current_AC = self.AC_data

    def apply_IDP(self, IDPs_to_filter):
        self.current_data= self.metadata
        self.current_AC = self.AC_data
        if (len(IDPs_to_filter)>0):
            self.current_data = self.current_data[self.current_data["IDP name"].isin(IDPs_to_filter)]
        target_genes = self.current_data["Target_Gene"]
        transcription_genes = self.current_data["Transcription_Factor"]
        self.current_AC = self.current_AC[(self.current_AC.iloc[:,0].isin (target_genes)) | (self.current_AC.iloc[:,0].isin (transcription_genes))]
        self.current_BS = self.current_BS[(self.current_BS.iloc[:,0].isin (target_genes)) | (self.current_BS.iloc[:,0].isin (transcription_genes))]


    def circos_setup(self):
        target_scatter=[]
        for i, row in self.current_data.iterrows():
            color="#66CCFF"
            promoter_text= "Enhancer"
            if row["Promoter"]==True:
                color="#9933FF"
                promoter_text= "Promoter"
            data = {'name':row["Transcription_Factor"],
            'fusion': promoter_text+"\n"+"Transcription Factor: "+row["Transcription_Factor"]+"\n"+"Target Gene: "+row["Target_Gene"]+"\n"+"Chromosome: "+str(row["chr"])+"\n"+"Start: "+str(row["Start"])+"\n"+"Stop: "+str(row["Stop"]),
             'g1chr': str(row["TGchr"]),'g1start': str(row["TGstart"]), 'g1end': str(row["TGend"]), 'g1name': str(row["Target_Gene"]),
            'g2chr': str(row["TFchr"]), 'g2start': str(row["TFstart"]), 'g2end': str(row["TFend"]), 'g2name': str(row["Transcription_Factor"]),
            'color':color}
            #data={'chr': row['chr'], 'nonref':row['nonref allele'],'pos':row['pos'],'ref':row['ref allele'],'start': row["Start"] , 'end': row["Stop"], 'name': row["Target_Gene"] }
            if (row["TGstart"]!=-1 and row["TGend"]!=-1 and row["TFstart"]!=-1 and row["TFend"]!=-1):
                target_scatter.append(data)

        return {'Points':target_scatter}

    def gsea_enrichement(self):
        matplotlib.use('Agg')
        gene_sets = ['GO_Biological_Process_2015','GO_Cellular_Component_2015','GO_Molecular_Function_2015','KEGG_2015']

        for gs in gene_sets:
            if os.path.exists("./static/images/"+gs+"_bar_enchr.png"):
              os.remove("./static/images/"+gs+"_bar_enchr.png")
            else:
              print("-------------------"+gs+"_bar_enchr.png does not exist----------------------------")

        genes_to_test = ((self.current_data["Transcription_Factor"]).append(self.current_data["Target_Gene"])).drop_duplicates()
        output = gseapy.enrichr(gene_list = genes_to_test, description='pathway', gene_sets=gene_sets, outdir='./static/images',cutoff=0.05, format='png')
        output.results["Label"]=output.results["Term"].str.split("(", n = 1, expand = True)[0]

        for g in output.results["Gene_set"].unique():
            gene_set=output.results[output.results["Gene_set"]==g]
            gene_set = gene_set.reset_index()
            gene_set["Label"]=gene_set["Term"].str.split("(", n = 1, expand = True)[0]

            gene_set = gene_set[gene_set["Genes"].str.split(";").str.len()>1]
            if not gene_set.empty:
                num_genes = gene_set["Genes"].str.split(";").str.len()
                p_values = 0-np.log(gene_set["Adjusted P-value"])

                cmap = matplotlib.cm.get_cmap('YlOrRd')
                norm = matplotlib.colors.Normalize(vmin=p_values.min(), vmax=p_values.max())

                colors=[]
                gene_set = gene_set.reset_index(drop=True)
                for i in p_values:
                    rgba = cmap(norm(i))
                    colors.append(rgba)
                plt.figure(figsize=(10,10))
                plt.title(g,fontsize=15)
                plt.xticks(fontsize=len(num_genes)/20)
                plt.xlabel("Number of Genes", fontsize=15)
                plt.yticks(gene_set.index,gene_set["Label"],fontsize=15)
                plt.barh(gene_set.index, num_genes, color=colors)
                # Add a colorbar

                sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
                sm.set_array([])

                cbar = plt.colorbar(sm)
                cbar.set_label('-log Adjusted P-value', rotation=270,labelpad=25, fontsize=15)

                plt.savefig("./static/images/"+g+"_bar_enchr.png",bbox_inches="tight")


    def AC_expression_setup(self):
        self.current_AC = self.current_AC.set_index('Gene')
        df = (self.current_AC).T
        print(df.mean(),df.max())
        df_norm_col=(df-df.mean())/(df.max() - df.min())
        plt.figure(figsize=(10,10))
        fig = sns.heatmap(df_norm_col, cmap='YlOrRd')
        plt.savefig("./static/images/ac_gene.png",bbox_inches="tight")

    def brainspan_setup(self):
        self.current_BS= self.current_BS.set_index('Gene')
        df = self.current_BS.T
        df_norm_col=(df-df.mean())/df.std()
        plt.figure(figsize=(10,10))
        fig = sns.heatmap(df_norm_col, cmap='YlOrRd')
        plt.savefig("./static/images/bs_gene.png",bbox_inches="tight")
