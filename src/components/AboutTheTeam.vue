<template>
  <div class="chart-container">
    <div ref="chart" class="chart"></div>
  </div>
</template>

<script setup>
import { ref, onMounted } from 'vue';
import * as d3 from 'd3';

const props = defineProps({
  data: {
    type: Array,
    required: true,
    validator(value) {
      return value.every(
        item => 'name' in item && 'link' in item && 'image' in item
      );
    }
  }
});

const chart = ref(null);

onMounted(() => {
  createChart();

  window.addEventListener('resize', createChart);

  function createChart() {
    const container = chart.value;
    const data = [...props.data]; //shuffleArray([...props.data]); // Randomize data order
    // Clear any existing SVG elements
    d3.select(container).selectAll('*').remove();

    const margin = { top: 20, right: 20, bottom: 20, left: 20 };
    const width = container.clientWidth - margin.left - margin.right;
    const radius = 100; // Fixed circle size
    const padding = 10; // Padding between circles

    // Calculate the number of columns and rows based on screen width
    const columns = Math.floor(width / (2 * radius));
    const rows = Math.ceil(data.length / columns);
    
    // Calculate the required height dynamically
    const height = rows * (2 * radius + padding) + margin.top + margin.bottom;

    const svg = d3
      .select(chart.value)
      .append('svg')
      .attr('width', width + margin.left + margin.right)
      .attr('height', height)
      .attr('viewBox', `0 0 ${width + margin.left + margin.right} ${height}`)
      .append('g')
      .attr('transform', `translate(${margin.left},${margin.top})`);

    // Build array to use to position elements based on indices
    let xCenter = [];
    for (let i = 0; i <= data.length - 1; i++) {
      xCenter.push(i * 100)
    }
    d3.forceSimulation(data)
      .force("x", d3.forceX((d, i) => {
          return xCenter[i];
      }).strength(1))
      .force("center", d3.forceCenter(width / 2, (height - margin.top - margin.bottom) / 2).strength(1))
      .force("charge", d3.forceManyBody().strength(30))

      .force("collision", d3.forceCollide().radius(radius + padding))
      .on("tick", ticked);

    const node = svg
      .selectAll("g")
      .data(data)
      .enter()
      .append("g")
      .attr('class', 'node')
      .on('click', function (event, d) {
        console.log('Clicked:', d.link); // Log to ensure click event is triggered
      });

    node
      .append('clipPath')
      .attr('id', (d, i) => `clip-${i}`)
      .append('circle')
      .attr('r', radius);

    node
      .append('circle')
      .attr('class', 'face')
      .attr('r', radius)
      .attr('fill', 'black')
      .attr('stroke-width', '5px')
      .attr('clip-path', (d, i) => `url(#clip-${i})`);

    node
      .append('image')
      .attr('xlink:href', d => d.image)
      .attr('width', radius * 2.1)
      .attr('height', radius * 2.1)
      .attr('x', -1.05 * radius)
      .attr('y', -radius)
      .attr('clip-path', (d, i) => `url(#clip-${i})`);

    node
      .append('circle')
      .attr('class', 'overlay')
      .attr('r', radius)
      .attr('fill', 'rgba(0, 0, 0, 0.5)')
      .style('opacity', 0);

    node
      .append('text')
      .attr('class', 'name-text')
      .attr('x', 0)
      .attr('y', 0)
      .attr('dx', '0em')
      .attr('dy', '0em')
      .attr('text-anchor', 'middle')
      .attr("dominant-baseline", "central")
      .attr("text-width", radius * 2)
      .style('opacity', 0)
      .text(d => d.name)
      .call(d => wrap(d, {shift: true}));

    // append link
    node    
      .append("a")
        .attr("href", d => d.link? d.link : null)
        .attr("target", "_blank")
        .append("svg:rect")
          .attr("y", -radius)
          .attr("x", -radius)
          .attr("height", radius * 2)
          .attr("width", radius * 2)
          .style("fill", "transparent")
          .attr('clip-path', (d, i) => `url(#clip-${i})`);

    node
      .on('mouseover', function () {
        d3.select(this).select('.overlay')
          .transition()
          .duration(200)
          .style('opacity', 1);
        d3.select(this).select('.name-text')
          .transition()
          .duration(200)
          .style('opacity', 1);
      })
      .on('mouseout', function () {
        d3.select(this).select('.overlay')
          .transition()
          .duration(200)
          .style('opacity', 0);
        d3.select(this).select('.name-text')
          .transition()
          .duration(200)
          .style('opacity', 0);
      });

    function ticked() {
      node.attr("transform", d => {
        d.x = Math.max(radius, Math.min(width - radius, d.x));
        if (rows == 1) {
          d.y = (height - margin.top - margin.bottom) / 2;
        } else {
          d.y = Math.max(radius, Math.min(height - margin.top - margin.bottom - radius, d.y));
        }
        return `translate(${d.x},${d.y})`;
      });
    }
  }

  // function shuffleArray(array) {
  //   for (let i = array.length - 1; i > 0; i--) {
  //     const j = Math.floor(Math.random() * (i + 1));
  //     [array[i], array[j]] = [array[j], array[i]];
  //   }
  //   return array;
  // }

  // https://gist.github.com/mbostock/7555321
  function wrap(text, {
      shift = false
  }) {
      text.each(function() {
          // split at hyphens, too: text.text().split(/\s|-+/) (but they don't get re-inserted)
          var text = d3.select(this),
          words = text.text().split(/\s+/).reverse(),
          word,
          line = [],
          lineNumber = 0,
          lineHeight = 1.1, // ems
          width = text.attr("text-width"),
          baseline = text.attr("dominant-baseline"),
          x = text.attr("x"),
          y = text.attr("y"),
          dy = parseFloat(text.attr("dy")),
          dx = parseFloat(text.attr("dx")),
          tspan = text.text(null).append("tspan").attr("y", y).attr("dy", dy + "em").attr("dominant-baseline", baseline);;

          while ((word = words.pop())) {
          line.push(word);
          tspan.text(line.join(" "));
              if (tspan.node().getComputedTextLength() > width) {
              line.pop();
              tspan.text(line.join(" "));
              line = [word];
              tspan = text.append("tspan").attr("x", x).attr("y", y).attr("dx", dx).attr("dy", ++lineNumber * lineHeight + dy + "em").attr("dominant-baseline", baseline).text(word);
              }
          }

          // https://stackoverflow.com/questions/60558291/wrapping-and-vertically-centering-text-using-d3-js
          if (lineNumber > 0  && shift) {
              const startDy = -(lineNumber * (lineHeight / 2));
              text
                  .selectAll("tspan")
                  .attr("dy", (d, i) => startDy + lineHeight * i + "em");
          }
      }
  )};
});
</script>


<style scoped lang="scss">
.chart-container {
  display: flex;
  justify-content: center;
  align-items: center;
  width: 100%;
}

.chart {
  width: 100%;
  height: 100%;
}
</style>
<style lang="scss">
/* css for elements added/classed w/ d3 */
.name-text {
  font-weight: 700;
  fill: #ffffff;
  pointer-events: 'none';
}
</style>