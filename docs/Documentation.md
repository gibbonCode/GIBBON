---
layout: pagefullwidth
title: "Documentation"
logo: "img/home-bg.jpg"
description: "Help and demo links"
header-img: "img/home-bg.jpg"
ordernumber: 4
---

# Documentation
__The documentation is a work in progress__. Not all functions have associated help files and not all functionality is covered by the demos. The help files currently cover about 50% of the functions and the demos mainly cover the use of FEBio.  

#### Table of content
* [The MATLAB integrated help](#helpMatlab)
* [Function help files](#help)
* [Demo files](#demo)

## The MATLAB integrated help <a name="helpMatlab"></a>  
Follow the installation instructions to integrate and access GIBBON documentation from within MATLAB. The name for all function help files (the files that generate the help/documentation when published with MATLAB) starts with `HELP_`, all demo files start with `DEMO_`. This way users may explore/open/edit these files by typing `open HELP_functionName` or `open DEMO_functionName` in the MATLAB command window

## Function help files <a name="help"></a>

<div>
  {% for help_html in site.static_files %}
     {% if help_html.path contains '/html' %}
        {% if help_html.basename contains 'HELP' %}
        {% if help_html.extname contains 'html' %}            

        <details>
          <summary><u><i class="fa fa-chevron-down" aria-hidden="true"></i> {{help_html.basename}}</u></summary>
          <span style="font-size:65%;">                    
            <iframe   src="{{ site.baseurl}}{{ help_html.path}}" width="100%" height="1000" frameborder="1" allowfullscreen> </iframe>
          </span>
        </details>
                  {% for image in site.static_files %}
                    {% if image.path contains help_html.basename %}
                       {% if image.path contains '.png' %}
                          {% if image.basename != help_html.basename %}
                            <img src="{{ site.baseurl }}{{ image.path }}" width="20%">  
                         {% endif %}              
                       {% endif %}              
                    {% endif %}              
                  {% endfor %}
                  <hr>  

        {% endif %}
        {% endif %}
     {% endif %}
 {% endfor %}

</div>

## Demo files <a name="demo"></a>  
<div>
{% for help_html in site.static_files %}
   {% if help_html.path contains '/html' %}
      {% if help_html.basename contains 'DEMO' %}
      {% if help_html.extname contains 'html' %}            

      <details>
        <summary><u><i class="fa fa-chevron-down" aria-hidden="true"></i> {{help_html.basename}}</u></summary>
        <span style="font-size:65%;">                    
          <iframe   src="{{ site.baseurl}}{{ help_html.path}}" width="100%" height="1000" frameborder="1" allowfullscreen> </iframe>
        </span>
      </details>
                {% for image in site.static_files %}
                  {% if image.path contains help_html.basename %}
                     {% if image.path contains '.png' %}
                        {% if image.basename != help_html.basename %}
                          <img src="{{ site.baseurl }}{{ image.path }}" width="20%">  
                       {% endif %}              
                     {% endif %}              
                  {% endif %}              
                {% endfor %}
                <hr>  

      {% endif %}
      {% endif %}
   {% endif %}
{% endfor %}
</div>
