---
layout: pagefullwidth
title: "Documentation"
logo: "img/home-bg.jpg"
description: "Help and demo links"
header-img: "img/home-bg.jpg"
ordernumber: 4
---

# Documentation
If GIBBON is installed properly the documentation is also available (and fully searchable) from within MATLAB (see [installation instructions]({{ site.baseurl }}/Installation/)). The name for all function help files (the files that generate the help/documentation when published with MATLAB) starts with `HELP_`, all demo files start with `DEMO_`. This way users may explore/open/edit these files by typing `open HELP_functionName` or `open DEMO_demoName` in the MATLAB command window.

__The documentation is a work in progress__. Not all functions have associated help files and not all functionality is covered by the demos. The help files currently cover about 50% of the functions and the demos mainly cover the use of FEBio.  

# Jump to sections:
# <i class="fa fa-arrow-right" aria-hidden="true"></i> [Function help files](#help)
# <i class="fa fa-arrow-right" aria-hidden="true"></i> [Demo files](#demo)

# Function help files <a name="help"></a>

<div style="background-color:lightgray">
  {% for help_html in site.static_files %}
    {% if help_html.path contains '/html' %}
      {% if help_html.basename contains 'HELP' %}
        {% if help_html.extname contains 'html' %}            

        <details>
          <summary>
            <u><i class="fa fa-plus-square" aria-hidden="true"></i> {{help_html.basename}}</u>
            <br>
            {% for image in site.static_files %}
              {% if image.path contains help_html.basename %}
                 {% if image.path contains '.jpg'%}
                    {% if image.basename != help_html.basename %}
                      <img src="{{ site.baseurl }}{{ image.path }}" width="5%">  
                    {% endif %}              
                 {% endif %}              
              {% endif %}              
            {% endfor %}
          </summary>

          <a href="{{ site.baseurl}}{{ help_html.path}}"> <i class="fa fa-arrow-right" aria-hidden="true"></i>  <i class="fa fa-file" aria-hidden="true"></i> <u>Full page view</u> </a>     

          <span style="font-size:65%;">            
            <iframe   src="{{ site.baseurl}}{{ help_html.path}}" width="100%" height="1000" frameborder="1" allowfullscreen> </iframe>
          </span>
        </details>

        <hr>

        {% endif %}
      {% endif %}
    {% endif %}
  {% endfor %}  
</div>

# Demo files <a name="demo"></a>  
<div style="background-color:lightgray">
  {% for help_html in site.static_files %}
    {% if help_html.path contains '/html' %}
      {% if help_html.basename contains 'DEMO' %}
        {% if help_html.extname contains 'html' %}            

        <details>
          <summary>
            <u><i class="fa fa-plus-square" aria-hidden="true"></i> {{help_html.basename}}</u>
            <br>
            {% for image in site.static_files %}
              {% if image.path contains help_html.basename %}
                 {% if image.path contains '.jpg'%}
                    {% if image.basename != help_html.basename %}
                      <img src="{{ site.baseurl }}{{ image.path }}" width="5%">  
                    {% endif %}              
                 {% endif %}              
              {% endif %}              
            {% endfor %}
          </summary>

          <a href="{{ site.baseurl}}{{ help_html.path}}"> <i class="fa fa-arrow-right" aria-hidden="true"></i>  <i class="fa fa-file" aria-hidden="true"></i> <u>Full page view</u> </a>     

          <span style="font-size:65%;">            
            <iframe   src="{{ site.baseurl}}{{ help_html.path}}" width="100%" height="1000" frameborder="1" allowfullscreen> </iframe>
          </span>
        </details>

        <hr>

        {% endif %}
      {% endif %}
    {% endif %}
  {% endfor %}  
</div>
