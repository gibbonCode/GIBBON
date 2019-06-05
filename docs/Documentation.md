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

## HELP and DEMO files
Click the links below to expand the documentation. Click on `Full page view` to open the documentation file in a seperate window.   

<div style="background-color:lightgray">

<style>         
   table {
         border-collapse: collapse;
         width: 100%;
         }

   td, th {
         border: 1px solid #dddddd;
         text-align: left;
         }

   tr:nth-child(even) {
         background-color: #dddddd;
         }
</style>

<table style="width:100%">
 <tr>
   <th> </th>
 </tr>

   {% for help_html in site.static_files %}
     {% if help_html.path contains '/html' %}
       {% if help_html.basename contains 'HELP'%}
         {% if help_html.extname contains 'html' %}            
         <tr>
           <td>
         <details>
           <summary>
             <u> <i class="fa fa-plus-square" aria-hidden="true"></i> {{help_html.basename}} </u>
             <br>
           </summary>
           <a href="{{ site.baseurl}}{{ help_html.path}}"> <i class="fa fa-arrow-right" aria-hidden="true"></i>  <i class="fa fa-file" aria-hidden="true"></i> <u>Full page view</u> </a>     
           <span style="font-size:65%;">            
             <iframe   src="{{ site.baseurl}}{{ help_html.path}}" width="100%" id="chgtext" height="2500" frameborder="1" scrolling="yes" style="font-size:65%;"> </iframe>
           </span>
         </details>
         </td>   
       </tr>
         {% endif %}
       {% endif %}
     {% endif %}
   {% endfor %}  

   {% for help_html in site.static_files %}
     {% if help_html.path contains '/html' %}
       {% if help_html.basename contains 'DEMO' %}
         {% if help_html.extname contains 'html' %}            
         <tr>
           <td>
         <details>
           <summary>
             <u> <i class="fa fa-plus-square" aria-hidden="true"></i> {{help_html.basename}} </u>
             <br>
           </summary>
           <a href="{{ site.baseurl}}{{ help_html.path}}"> <i class="fa fa-arrow-right" aria-hidden="true"></i>  <i class="fa fa-file" aria-hidden="true"></i> <u>Full page view</u> </a>     
           <span style="font-size:65%;">            
             <iframe   src="{{ site.baseurl}}{{ help_html.path}}" width="100%" id="chgtext" height="2500" frameborder="1" scrolling="yes" style="font-size:65%;"> </iframe>
           </span>
         </details>
         </td>   
       </tr>
         {% endif %}
       {% endif %}
     {% endif %}
   {% endfor %}  
</table>

</div>
