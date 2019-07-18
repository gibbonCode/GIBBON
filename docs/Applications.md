---
layout: pagefullwidth
title: "Applications"
description: ""
header-img: "img/home-bg.jpg"
ordernumber: 9
---

<style>

.card {
    box-shadow: 0 8px 16px 0 rgba(0,0,0,0.2);
    transition: 0.3s;
    width: 100%;
		height: 350px;
    overflow: hidden;		
    border-style: solid;
    border-color: gray;
    border-width: thin;
    border-radius: 10px; /* 5px rounded corners */
    background-color: white;
    margin-bottom: 30px;
}
.card-text {
    word-wrap: break-word;
    margin-left: 3rem;
}

.card:hover {
    box-shadow: 0 16px 32px 0 rgba(0,0,255,0.2);
    z-index: 2;
-webkit-transition: all 200ms ease-in;
-webkit-transform: scale(1.05);
-ms-transition: all 200ms ease-in;
-ms-transform: scale(1.05);
-moz-transition: all 200ms ease-in;
-moz-transform: scale(1.05);
transition: all 200ms ease-in;
transform: scale(1.05);
}


</style>

# Projects using GIBBON

<div class="row">

{% for application in site.applications %}
<div class="col-xs-12 col-sm-6 col-md-3 col-l-3 col-xl-3 d-flex align-items-stretch">
  <a href="{{ application.url | prepend: site.baseurl }}" >
	<div class="card" align="center">

					<b class="post-title">{{ application.title }}</b>

					{% if application.subtitle %}
        	<br />
	          <i class="post-subtitle">{{ application.subtitle }}</i>
        	{% endif %}

					{% if application.authors %}					
					<br />
							<i>{{ application.authors }}</i>						
					{% endif %}

					{% if application.thumbnail %}
        	<br />
        	    <img alt="" src="{{ site.baseurl }}/img/applications/{{ application.thumbnail }}" style="width:100%; max-height: 220px; max-width: 220px">
        	{% endif %}

	</div>
  </a>
</div>
{% endfor %}

</div>
