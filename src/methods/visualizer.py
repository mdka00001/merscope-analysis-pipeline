from src.methods.dependencies import *
from src.methods.base import BaseMethod
from observable_jupyter import embed

class Visualizer(BaseMethod):
    def run(self):
        """
        Perform clustering on the input data using Louvain and Leiden algorithms.
        """
        # Load the input data
        print(f"Loading data from {self.adata}")
        adata = sc.read_h5ad(self.adata)

        # initialize meta_gene
        cell_by_gene = pd.read_csv(self.input_cell_by_gene, index_col=0)
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

        
        # hierarchically cluster leiden signatures
        net = Network(CGM2)
        net.load_df(sig_leiden_clip, meta_col=meta_leiden, col_cats=['leiden', 'Cell_Type', 'cell counts'], 
                    meta_row=meta_gene, row_cats=['Markers'])
        net.filter_threshold(0.01, axis='row')
        net.set_global_cat_colors(df_colors)
        net.cluster()

        # Collect data for interactive visualization
        ################################################
        cell_by_gene.index = range(len(cell_by_gene.index.tolist()))
        gex_int = cell_by_gene.astype(int)

        # new gex_dict
        ###################
        gex_dict = {}
        min_exp = 5
        max_exp = 255
        max_pos_cells = 40000
        for inst_gene in gex_int.columns.tolist():
            if 'Blank' not in inst_gene:
                ser_gene = gex_int[inst_gene]
                ser_gene[ser_gene > max_exp] = max_exp
                ser_gene = ser_gene[ser_gene > min_exp]
                ser_gene = ser_gene.astype(np.uint8)    

                if ser_gene.shape[0] > max_pos_cells:
                # print(inst_gene, ser_gene.shape)
                    keep_cells = ser_gene.sort_values(ascending=False).index.tolist()[:max_pos_cells]
                    ser_gene = ser_gene[keep_cells]

                # dictonary
                gex_dict[inst_gene] = ser_gene.to_dict()

        zip_gex_dict = json_zip(gex_dict)

        
        df_pos = adata.obs[['center_x', 'center_y', 'leiden']]
        df_pos[['center_x', 'center_y']] = df_pos[['center_x', 'center_y']].round(2)
        df_pos.columns = ['x', 'y', 'leiden']
        df_pos['y'] = -df_pos['y']
        df_umap = adata.obsm.to_df()[['X_umap1', 'X_umap2']].round(2)
        df_umap.columns = ['umap-x', 'umap-y']

        df_name = pd.DataFrame(df_pos.index.tolist(), index=df_pos.index.tolist(), columns=['name'])

        df_obs = pd.concat([df_name, df_pos, df_umap], axis=1)
        data = df_obs.to_dict('records')

        obs_data = {
            # 'gex_dict': gex_dict,
            # 'gex_dict': {},
            'data': data, 
            'cat_colors': cat_colors,
            'network': net.viz    
        }

        zip_obs_data = json_zip(obs_data)


        inputs = {
            'zoom': -3.75, 
            'ini_cat': 'leiden',
            'ini_map_type': 'Spatial',
            'ini_min_radius': 1.0,
            'zip_obs_data': zip_obs_data,
            'zip_gex_dict': zip_gex_dict,
            'gex_opacity_contrast_scale': 0.85, 
            'center_x': 5000,
            'center_y': 5000    
        }

        # initialize df_points DataFrame for neighborhood calculations
        df_points = pd.concat([df_name, df_pos], axis=1)

        embed('@vizgen/umap-spatial-heatmap-single-cell-0-3-1', cells=['viewof cgm', 'dashboard'], inputs=inputs, display_logo=False)



        