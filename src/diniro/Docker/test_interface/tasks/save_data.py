import scanpy as sc
import random
import pandas as pd
import json
from users.models import UserFiles


class SaveData:

    def read_data(userrr, map='non', obs='non'):
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        print('START task save_data.py')

        if map == 'non' and obs == 'non':
            current_user = UserFiles.objects.get(user=userrr)
            #(current_user.main_anndata)
            print(current_user.main_anndata)
            anndata = sc.read_h5ad(current_user.main_anndata)
            plot_maps = list(anndata.obsm)
            # if 'X_umap' in plot_maps:
            #     x, y = anndata.obsm['X_umap'][:, 0], anndata.obsm['X_umap'][:, 1]
            # else:
            x, y = anndata.obsm[plot_maps[0]][:, 0], anndata.obsm[plot_maps[0]][:, 1]
            observations = list(anndata.obs)
            # if 'louvain' in observations:
            #     c = anndata.obs['louvain']
            # else:
            c = anndata.obs[observations[0]]
            c = list(c)

        elif map != 'non' and obs != 'non':

            current_user = UserFiles.objects.get(user=userrr)
            anndata = sc.read_h5ad(current_user.main_anndata)
            plot_maps = list(anndata.obsm)
            if map == 'X_diffmap':
                x, y = anndata.obsm[str(map)][:, 1], anndata.obsm[str(map)][:, 2]
            else:
                x, y = anndata.obsm[str(map)][:, 0], anndata.obsm[str(map)][:, 1]
            observations = list(anndata.obs)
            c = anndata.obs[obs]
            c = list(c)
            #print(anndata.obs.groupby([obs]).apply(len))

        N_Cells = len(list(anndata.obs_names))
        N_Genes = len(list(anndata.var_names))
        N_Groups = len(list(anndata.obsm))

        original_colors = c

        def cat_to_color(categorical):
            new_colors = []
            colors = ["#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                        "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                        "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                        "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                        "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
                        "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                        "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                        "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                        "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                        "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
                        "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                        "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                        "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C"]
            size = len(set(categorical))
            colo = random.sample(colors, k=size)
            for x in categorical:
                new_colors.append(colo[x])
            return new_colors

        from sklearn.preprocessing import LabelEncoder
        label_encoder = LabelEncoder()
        c_categorical = label_encoder.fit_transform(c)

        df = pd.DataFrame({'x': x, 'y': y})
        result = df.to_json(orient='records')
        points = json.loads(result)

        if len(list(set(c))) < 64:
            cc = cat_to_color(c_categorical)
            color = cc
        else:
            color = c

        model, created = UserFiles.objects.get_or_create(user=userrr)
        model.points = points
        model.colors = color
        model.save()

        print("DONE task save_data.py")
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        return N_Cells, N_Genes, N_Groups, points, c, color, plot_maps, observations

    def legend_check(unique_leb, points, leb, cc):
        plotting_data = []
        if len(unique_leb) < 64:

            for val in unique_leb:
                res_index = [i for i, value in enumerate(leb) if value == val]
                thislist = []
                thislist1 = []
                thislist2 = []
                for index in res_index:
                    thislist.append(points[index])
                    thislist1.append(leb[index])
                    thislist2.append(cc[index])
                x = []
                y = []
                for v in thislist:
                    x.append([(k, v) for k, v in v.items()][0][1])
                    y.append([(k, v) for k, v in v.items()][1][1])
                plotting_data.append([x, y, thislist1[0], thislist2[0]])
        else:
            x = []
            y = []
            for v in points:
                x.append([(k, v) for k, v in v.items()][0][1])
                y.append([(k, v) for k, v in v.items()][1][1])
            plotting_data.append([x, y, 'data', cc])

        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        return plotting_data
