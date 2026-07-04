---
layout: pagefullwidth
title: "Documentation"
logo: "img/home-bg.jpg"
description: "Help and demo links"
header-img: "img/home-bg.jpg"
ordernumber: 4
---

# Documentation
If GIBBON is installed properly the documentation is also available (and fully searchable) from within MATLAB (see [installation instructions]({{ site.baseurl }}/Installation/)).

__The documentation is a work in progress__. Not all functions have associated help files and not all functionality is covered by the demos. The help files currently cover about 50% of the functions and the demos mainly cover the use of FEBio.  

## HELP and DEMO files
The name for all function help files starts with `HELP_`, and all demo files start with `DEMO_`. This way users may explore/open/edit these files by typing `open HELP_functionName` or `open DEMO_demoName` in the MATLAB command window.

<div style="background-color:white">

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

         {
           box-sizing: border-box;
         }

         #myInput {
           background-image: url('/img/icons/search.png');
           background-position: 10px 12px;
           background-repeat: no-repeat;
           width: 100%;
           font-size: 16px;
           padding: 12px 20px 12px 40px;
           border: 1px solid #ddd;
           margin-bottom: 12px;
         }

         #myUL {
           list-style-type: none;
           padding: 0;
           margin: 0;
         }

         #myUL li a {
           border: 1px solid #ddd;
           margin-top: -1px; /* Prevent double borders */
           background-color: #f6f6f6;
           padding: 12px;
           text-decoration: none;
           font-size: 18px;
           color: black;
           display: block
         }

         #myUL li a:hover:not(.header) {
           background-color: #eee;
         }

</style>

<script>
function myFunction() {
    var input, filter, ul, li, a, i, txtValue;
    input = document.getElementById("myInput");
    filter = input.value.toUpperCase();
    ul = document.getElementById("myUL");
    li = ul.getElementsByTagName("li");
    for (i = 0; i < li.length; i++) {
        a = li[i].getElementsByTagName("a")[0];
        txtValue = a.textContent || a.innerText;
        if (txtValue.toUpperCase().indexOf(filter) > -1) {
            li[i].style.display = "";
        } else {
            li[i].style.display = "none";
        }
    }
}
</script>

<input type="text" id="myInput" onkeyup="myFunction()" placeholder=" Search documentation files.." title="Type in a name">


<ul id="myUL">

{% assign sortedFiles = site.static_files | sort_natural: 'basename' %}

   {% for help_html in sortedFiles %}
     {% if help_html.path contains '/html' %}
       {% if help_html.basename contains 'HELP'%}
         {% if help_html.extname contains 'html' %}            
            <li> <a href="{{ site.baseurl}}{{ help_html.path}}"> {{help_html.basename  | remove: "HELP_"}} </a>  </li>        
         {% endif %}
       {% endif %}
     {% endif %}
   {% endfor %}  

  {% for help_html in sortedFiles %}
     {% if help_html.path contains '/html' %}
       {% if help_html.basename contains 'DEMO' %}
         {% if help_html.extname contains 'html' %}            
            <li> <a href="{{ site.baseurl}}{{ help_html.path}}"> {{help_html.basename}} </a>  </li>        
         {% endif %}
       {% endif %}
     {% endif %}
   {% endfor %}  

</ul>

</div>
