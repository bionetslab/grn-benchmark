# def about_remove(request):
#     print('************Here*************')
#     current_user = UserFiles.objects.get(user=request.user.username)
#     anndata = sc.read_h5ad(current_user.main_anndata)
#     x, y = anndata.obsm['X_umap'][:, 0], anndata.obsm['X_umap'][:, 1]
#     c = anndata.obs['louvain']
#     c = list(c)
#
#     for i in range(len(c)):
#         c[i] = c[i].decode("utf-8")
#
#     def cat_to_color(categorical):
#         new_colors = []
#         colors = ['#FF0000', '#FFFF00', '#008000', '#0000FF', '#FF00FF', '#800080', '#FF00FF', '#00FF00', '#008080',
#                   '#00FFFF', '#0000FF']
#         size = len(set(categorical))
#         #print(size)
#         colo = random.choices(colors, k=size)
#         #print(colo)
#         for x in categorical:
#             new_colors.append(colors[x])
#         return new_colors
#
#     from sklearn.preprocessing import LabelEncoder
#     label_encoder = LabelEncoder()
#     c_categorical = label_encoder.fit_transform(c)
#     #print(c_categorical)
#
#     cc = cat_to_color(c_categorical)
#
#     df = pd.DataFrame({'x': x, 'y': y})
#     result = df.to_json(orient='records')
#     points = json.loads(result)
#
#     model, created = UserFiles.objects.get_or_create(user=request.user.username)
#     model.points = points
#     model.colors = cc
#     model.save()
#     contex = {'points': points,
#               'colors': cc}
#     return print("DONE")
#     #return render(request, 'test_interface/about.html', contex)


# < script >
#
# $(document).on('submit', '#post-form', function(e)
# {
# $.ajax({
#     type: 'POST',
#     url: '{% url "test-upload" %}',
#     data: {
#         map: $('#map_selection').val(),
#               obs: $('#obs_selection').val(),
#                     csrfmiddlewaretoken:$('input[name=csrfmiddlewaretoken]').val(),
#                                          action: 'post'
# },
# success: function(json)
# {
#     document.getElementById("post-form").reset();
#
# },
# });
# });
# < / script >

#< !button id = "sel"class ="btn btn-primary btn-lg" role="button" onclick="go2()" > Select < ! / button >
#    <!button type="submit" class="btn btn-primary" >Select<!/button>


# var
# fruits = [];
# fruits.push(x);
# document.getElementById("adds1").innerHTML = fruits.length;


# < script
# type = "text/javascript"
# src = "https://unpkg.com/vis-network/standalone/umd/vis-network.min.js"
# > < / script >
#
# < style
# type = "text/css" >
# .mynetwork
# {
#     width: 600px;
# height: 400
# px;
# border: 1
# px
# solid
# lightgray;
# }
# < / style >


# { %
# for i in n_plots %}
#
# < script
# type = "text/javascript" >
# // create
# an
# array
# with nodes
#     var
#     nodes = new
#     vis.DataSet({{nodes | index: i | safe}});
#
# // create
# an
# array
# with edges
#     var
#     edges = new
#     vis.DataSet({{edges | index: i | safe}});
#
# // create
# a
# network
# var
# container = document.getElementById("mynetwork{{i}}");
# var
# data = {
#     nodes: nodes,
#     edges: edges,
# };
# var
# options = {};
# var
# network = new
# vis.Network(container, data, options);
# < / script >
#
# { % endfor %}