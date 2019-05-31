---
layout: page
title: "People"
subheadline: ""
teaser: ""
permalink: "/people/"
header:
---

{% for author in site.data.authors %}
  <div class="image">
    {% if author.url %}
      <a target="_blank" href="{{ author.url }}">
    {% endif %}
    {% if author.image %}
      <img src="{{ author.image }}" alt="Homepage" height="100px">
    {% else %}
      <img src="/img/anonymous.png" alt="Homepage" height="100px">
    {% endif %}
    {% if author.url %}
      </a>
    {% endif %}
    <div class="desc">{{ author.name }}</div>
  </div>
{% endfor %}
