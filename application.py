from flask import Flask, render_template,request
import json
import model
import os
import glob

application = Flask(__name__)

global AC

@application.route('/',methods=['GET', 'POST'])
def show_init():
    global IDP_options
    global data_model
    data_model= model.Model("data/imaging_network.csv","data/UCLA.csv", "data/brainspan.csv")
    IDP_options = data_model.metadata["IDP name"].unique()
    return render_template('index.html', IDP_options=IDP_options)

@application.route('/init', methods=['POST'])
def init():
    data_model.current_data = data_model.metadata
    return ""

@application.route('/apply_filters', methods=['POST'])
def apply_filters():
    IDP_selected = json.loads(request.get_data())
    data_model.apply_IDP(IDP_selected)
    circos_dict = data_model.circos_setup()
    enriched = data_model.gsea_enrichement()
    data_model.AC_expression_setup()
    data_model.brainspan_setup()
    return json.dumps({"circos_dict":circos_dict,"current_table":data_model.current_data.to_html()})

@application.route('/reset', methods = ['POST'])
def reset():
    data_model.reset_filter()
    return ""

@application.route('/descr', methods = ['POST'])
def read_files():
    IDP_selected = json.loads(request.get_data())["selected"]
    keywords = set(IDP_selected.split())
    descr_dict={}
    for word in keywords:
        file_path = 'data/descr/'+word+".txt"
        exists = os.path.isfile(file_path)
        if exists:
            file = open(file_path,'r')
            text = file.read()
            descr_dict[word] = [text]
            file.close()
            imges_file = (glob.glob('static/images/desc_images/'+word+"*"))
            if len(imges_file)>0:
                img_file= imges_file[0]
                descr_dict[word].append(img_file)
    return json.dumps(descr_dict)


if __name__ == '__main__':
    application.run()
