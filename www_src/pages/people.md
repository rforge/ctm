---
layout: page
title: "People"
permalink: "people/"
---

<div class="row">
{% for author in site.data.authors %}
  <div class="small-3 columns">
    {% if author.url %}
      <a target="_blank" href="{{ author.url }}">
    {% endif %}
    {% if author.image %}
    <img src="{{ site.urlimg }}{{ author.image }}" alt="Homepage" style="height:120px">
    {% else %}
    <img src="{{ site.urlimg }}anonymous.png" alt="Homepage"  style="height:120px">
    {% endif %}
    {% if author.url %}
      </a>
    {% endif %}
    <div class="desc">{{ author.name }}</div>
  </div>
{% endfor %}
</div>
