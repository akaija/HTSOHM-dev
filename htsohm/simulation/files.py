import os
from string import Template

import htsohm

def load_and_subs_template(template_name, params):
    htsohm_dir = os.path.dirname(htsohm.__file__)
    input_path = os.path.join(htsohm_dir, 'simulation', template_name)
    with open(input_path) as input_file:
        template = Template(input_file.read())
    return template.substitute(params)
