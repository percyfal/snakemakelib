# Copyright (C) 2015 by Per Unneberg
import os
from snakemakelib.config import sml_base_path
from bokeh.resources import CDN
from bokeh.templates import RESOURCES
from bokeh.embed import components
from bokeh.models.widget import Widget
from bokeh.util.string import encode_utf8

def static_html(template, resources=CDN, as_utf8=True, **kw):
    """Render static html document.

    This is a minor modification of bokeh.embed.file_html

    Args:
      template (Template): jinja2 HTML document template
      resources (Resources): a resource configuration for BokehJS assets
      kw: keyword argument list of bokeh components. Keywords must match with keywords in template
      
    Returns:
      html : standalone HTML document with embedded plot
    """
    plot_resources = RESOURCES.render(
        js_raw = resources.js_raw,
        css_raw = resources.css_raw,
        js_files = resources.js_files,
        css_files = resources.css_files,
    )

    with open (os.path.join(sml_base_path(), 'static/basic.css')) as fh:
        css = "".join(fh.readlines())

    # Hack to get on-the-fly double mapping. Can this be rewritten with e.g. map?
    def _update(kw):
        tmp = {}
        for k,v in kw.items():
            if (isinstance(v, Widget)):
                tmp.update({k : [{'script' : s, 'div' : d } for s,d in [components(v, resources)]][0]})
            elif (isinstance(v, dict)):
                if not v:
                    tmp.update(v)
                else:
                    v.update(_update(v))
            else:
                tmp.update({k:v})
        return tmp

    kw.update(_update(kw))
    if as_utf8:
        return encode_utf8(template.render(plot_resources=plot_resources, css=css, **kw))
    else:
        return template.render(plot_resources=plot_resources, css=css, **kw)
