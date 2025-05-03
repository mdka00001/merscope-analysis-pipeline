from src.methods.dependencies import *
from src.methods.base import BaseMethod
from src.output.outputDTO import OutputDTO

class ClusterAnnotation(BaseMethod):
    def run(self):
        """
        Perform clustering on the input data using Louvain and Leiden algorithms.
        """
        # Load the input data
        print(f"Loading data from {self.adata}")
        adata = sc.read_h5ad(self.adata)

        # initialize meta_gene
        meta_gene = pd.DataFrame(index=cell_by_gene.columns.tolist())

        # drop blanks for single cell analysis
        keep_genes = [x for x in cell_by_gene.columns.tolist() if 'Blank' not in x]

        cell_by_gene = cell_by_gene[keep_genes]
        meta_gene = meta_gene.loc[keep_genes]

        meta_gene['expression'] = cell_by_gene.sum(axis=0)

        vizgen_genes = meta_gene.index.tolist()

        ser_counts = adata.obs['leiden'].value_counts()
        ser_counts.name = 'cell counts'
        meta_leiden = pd.DataFrame(ser_counts)

        cat_name = 'leiden'
        sig_leiden = pd.DataFrame(columns=adata.var_names, index=adata.obs[cat_name].cat.categories)                                                                                                 
        for clust in adata.obs[cat_name].cat.categories: 
            sig_leiden.loc[clust] = adata[adata.obs[cat_name].isin([clust]),:].X.mean(0)
        sig_leiden = sig_leiden.transpose()
        leiden_clusters = ['Leiden-' + str(x) for x in sig_leiden.columns.tolist()]
        sig_leiden.columns = leiden_clusters
        meta_leiden.index = sig_leiden.columns.tolist()
        meta_leiden['leiden'] = pd.Series(meta_leiden.index.tolist(), index=meta_leiden.index.tolist())

        # generate colors for categories by plotting
        sc.pl.umap(adata, color="leiden", legend_loc='on data')
        cats = adata.obs['leiden'].cat.categories.tolist()
        colors = list(adata.uns['leiden_colors'])
        cat_colors = dict(zip(cats, colors))  

        # colors for clustergrammer2
        ser_color = pd.Series(cat_colors)
        ser_color.name = 'color'
        df_colors = pd.DataFrame(ser_color)
        df_colors.index = ['Leiden-' + str(x) for x in df_colors.index.tolist()]

        meta_gene['info'] = pd.Series('', index=meta_gene.index.tolist())
        df_colors.loc[''] = 'white'

        df_ref_panel_ini = pd.read_excel(self.ref_marker_panel, index_col=0)
        df_ref_panel = df_ref_panel_ini.iloc[1:,:1]
        df_ref_panel.index.name = None
        df_ref_panel.columns = ['Function']

        marker_genes = df_ref_panel[df_ref_panel['Function'].str.contains('marker')].index.tolist()
        common_marker_genes = list(set(meta_gene.index.tolist()).intersection(marker_genes))
        meta_gene.loc[common_marker_genes, 'Markers'] = df_ref_panel.loc[common_marker_genes, 'Function']
        meta_gene['Markers'] = meta_gene['Markers'].apply(lambda x: 'N.A.' if 'marker' not in str(x) else x)
        meta_gene['Markers'].value_counts()

        # Clip Z-score values for visual purposes
        sig_leiden_clip = deepcopy(sig_leiden)
        sig_leiden_clip[sig_leiden_clip < -5] = -5
        sig_leiden_clip[sig_leiden_clip >= 5] = 5

        # alphabetize genes
        sig_leiden_clip = sig_leiden_clip.loc[sorted(sig_leiden_clip.index.tolist())]
        sig_leiden_clip.head()

        meta_gene = pd.DataFrame(index=sig_leiden.index.tolist())
        meta_gene['Markers'] = pd.Series('N.A.', index=sig_leiden.index.tolist())
        meta_gene.loc[common_marker_genes, 'Markers'] = df_ref_panel.loc[common_marker_genes, 'Function']

        df_colors.loc['N.A.', 'color'] = 'white'

        meta_leiden['Cell_Type'] = pd.Series('N.A.', index=meta_leiden.index.tolist())
        num_top_genes = 30
        for inst_cluster in sig_leiden.columns.tolist():
            top_genes = sig_leiden[inst_cluster].sort_values(ascending=False).index.tolist()[:num_top_genes]
            
            inst_ser = meta_gene.loc[top_genes, 'Markers']
            inst_ser = inst_ser[inst_ser != 'N.A.']
            ser_counts = inst_ser.value_counts()

            max_count = ser_counts.max()

            max_cat = '_'.join(sorted(ser_counts[ser_counts == max_count].index.tolist()))
            max_cat = max_cat.replace(' marker', '').replace(' ', '-')

            print(inst_cluster, max_cat)  
            meta_leiden.loc[inst_cluster, 'Cell_Type'] = max_cat

        output = OutputDTO()
        output.add_data(adata)

        # Plotting
        plots = ["umap_leiden", "tSNE_leiden", "spatial_leiden"]

        for plot in plots:
            if plot == "tSNE_leiden":
                sc.pl.tsne(adata, color=["leiden"], show=False)
                fig = plt.gcf()
                output.add_plot("tSNE_leiden", fig)


            elif plot == "umap_leiden":
                sc.pl.umap(adata, color=["leiden"], show=False)
                fig = plt.gcf()
                output.add_plot("umap_leiden", fig)

            elif plot == "spatial_leiden":
                sc.pl.spatial(adata, color=["leiden"], show=False)
                fig = plt.gcf()
                output.add_plot("spatial_leiden", fig)

        return output