---
layout: page
title: "People"
subheadline: ""
teaser: ""
permalink: "/people/"
header:
---

{% for author in site.data.authors %}
  <div class="images">
    {% if author.url %}
      <a target="_blank" href="{{ author.url }}">
    {% endif %}
    {% if author.image %}
    <img src="../assets/{{ author.image }}" alt="Homepage" style="height:100px">
    {% else %}
      <img src="../assets/img/anonymous.png" alt="Homepage" style="height=100px">
    {% endif %}
    {% if author.url %}
      </a>
    {% endif %}
    <div class="desc">{{ author.name }}</div>
  </div>
{% endfor %}
