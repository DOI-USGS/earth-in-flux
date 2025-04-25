<template>
  <!---VizSection-->
  <VizSection :figures="true" :fig-caption="false">
    <!-- HEADING -->
    <template #heading> </template>
    <!-- FIGURES -->
    <template #aboveExplanation>
      <p v-html="text.paragraph1" />
      <FamilyInfoBox :activeFamily="activeFamily" />
    </template>

    <template #figures>
      <div class="chart-container single" ref="chart"></div>
    </template>
    <!-- FIGURE CAPTION -->
    <template #figureCaption> </template>
    <!-- EXPLANATION -->
    <template #belowExplanation>
      <p v-html="text.paragraph2" />
    </template>
  </VizSection>
</template>

<script setup>
import { onMounted, ref } from 'vue'
import * as d3 from 'd3'
import VizSection from '@/components/VizSection.vue'
import FamilyInfoBox from '@/components/FamilyInfoBox.vue'

// define props
const props = defineProps({
  text: { type: Object }
})

// global variables
const publicPath = import.meta.env.BASE_URL
const chart = ref(null)
const bodyCSS = window.getComputedStyle(document.body)
const lightBlue = bodyCSS.getPropertyValue('--light-blue')
const darkBlue = bodyCSS.getPropertyValue('--dark-blue')
const familyInfo = props.text.familyInfo
const defaultFamily = props.text.defaultFamily
const activeFamily = ref(props.text.defaultFamily)

// Declare behavior on mounted
// functions called here
onMounted(async () => {
  // Load the json data
  const data = await d3.json(publicPath + 'total_price.json')

  const processedData = lumpLowValueSpecies(data)

  // build chart
  buildChart(processedData)
})

function lumpLowValueSpecies(data, threshold = 500000) {
  const transformed = JSON.parse(JSON.stringify(data))

  transformed.children = transformed.children.map((family) => {
    const newChildren = []
    const countryMap = new Map()
    let lumpedSpeciesCount = 0

    family.children.forEach((species) => {
      const total = species.children.reduce((sum, country) => sum + country.value, 0)

      if (total < threshold) {
        lumpedSpeciesCount += 1
        species.children.forEach((country) => {
          if (countryMap.has(country.name)) {
            countryMap.set(country.name, countryMap.get(country.name) + country.value)
          } else {
            countryMap.set(country.name, country.value)
          }
        })
      } else {
        newChildren.push(species)
      }
    })

    if (countryMap.size > 0) {
      if (lumpedSpeciesCount === 1) {
        // Only one low-value species, push it through as is
        const onlySpecies = family.children.find((species) => {
          const total = species.children.reduce((sum, country) => sum + country.value, 0)
          return total < threshold
        })
        newChildren.push(onlySpecies)
      } else {
        // More than one , lump into "Other"
        const otherSpecies = {
          name: 'Other',
          children: [],
          lumpedSpeciesCount
        }

        countryMap.forEach((value, name) => {
          otherSpecies.children.push({ name, value })
        })

        newChildren.push(otherSpecies)
      }
    }

    return {
      ...family,
      children: newChildren,
      originalSpeciesCount: family.children.length,
      lumpedSpeciesCount
    }
  })

  return transformed
}

// https://observablehq.com/@d3/zoomable-circle-packing
function buildChart(data) {
  // Specify the chart’s dimensions.
  const width = 928
  const height = width

  // Create the color scale.
  const color = d3
    .scaleLinear()
    .domain([0, 5])
    .range([lightBlue, darkBlue])
    .interpolate(d3.interpolateHcl)

  // Compute the layout.
  const pack = (data) =>
    d3.pack().size([width, height]).padding(3)(
      d3
        .hierarchy(data)
        .sum((d) => d.value)
        .sort((a, b) => b.value - a.value)
    )
  const root = pack(data)

  // Create the SVG container.
  const svg = d3
    .select(chart.value)
    .append('svg')
    .attr('id', 'circle-pack-svg')
    .attr('viewBox', `-${width / 2} -${height / 2} ${width} ${height}`)
    .attr('width', width)
    .attr('height', height)

  // Append the nodes.
  const node = svg
    .append('g')
    .selectAll('circle')
    .data(root.descendants().slice(1))
    .join('circle')
    .attr('fill', (d) => (d.children ? color(d.depth) : 'white'))
    .attr('pointer-events', (d) => (!d.children ? 'none' : null))
    .on('mouseover', function () {
      d3.select(this).attr('stroke', '#000')
    })
    .on('mouseout', function () {
      d3.select(this).attr('stroke', null)
    })
    .on('click', (event, d) => {
      if (focus !== d) {
        zoom(event, d)
        event.stopPropagation()
      }
    })

  // Append the text labels.
  const label = svg
    .append('g')
    .attr('class', 'circle-packing-text')
    .attr('pointer-events', 'none')
    .attr('text-anchor', 'middle')
    .selectAll('text')
    .data(root.descendants())
    .join('text')
    .attr('class', 'fish-title')
    .style('fill-opacity', (d) => (d.parent === root ? 1 : 0))
    .style('stroke-opacity', (d) => (d.parent === root ? 1 : 0))
    .style('display', (d) => (d.parent === root ? 'inline' : 'none'))
    .text((d) => d.data.name)

  label
    .append('tspan')
    .attr('class', 'fish-value')
    .attr('dy', '0rem')
    .attr('dx', '0.3rem') // drop if dy set to 1rem
    // .attr("x", 0)    (if set dy to 1rem)
    .text((d) => d3.format('$.1s')(d.value).replace('G', 'B'))

  // Create the zoom behavior and zoom immediately in to the initial focus node.
  svg.on('click', (event, d) => {
    if (focus !== d) {
      zoom(event, d)
      event.stopPropagation()
    }
  })

  //  zoom-out on background click
  d3.select(chart.value)
    .select('svg')
    .on('click', (event) => {
      activeFamily.value = defaultFamily
      zoom(event, root)
    })

  let focus = root
  let view
  zoomTo([focus.x, focus.y, focus.r * 2.2])

  function zoomTo(v) {
    const k = width / v[2]

    view = v

    label.attr('transform', (d) => `translate(${(d.x - v[0]) * k},${(d.y - v[1]) * k})`)
    node.attr('transform', (d) => `translate(${(d.x - v[0]) * k},${(d.y - v[1]) * k})`)
    node.attr('r', (d) => d.r * k)
  }

  function zoom(event, d) {
    focus = d

    const transition = svg
      .transition()
      .duration(event.altKey ? 7500 : 750)
      .tween('zoom', () => {
        const i = d3.interpolateZoom(view, [focus.x, focus.y, focus.r * 2.2])
        return (t) => zoomTo(i(t))
      })

    label
      .filter(function (d) {
        return d.parent === focus || this.style.display === 'inline'
      })
      .transition(transition)
      .style('fill-opacity', (d) => (d.parent === focus ? 1 : 0))
      .style('stroke-opacity', (d) => (d.parent === focus ? 1 : 0))
      .on('start', function (d) {
        if (d.parent === focus) this.style.display = 'inline'
      })
      .on('end', function (d) {
        if (d.parent !== focus) this.style.display = 'none'
      })
    // format species level label
    function formatEconomicValue(val) {
      const suffixes = [
        { value: 1e9, symbol: 'B' },
        { value: 1e6, symbol: 'M' }
      ]

      for (let i = 0; i < suffixes.length; i++) {
        if (val >= suffixes[i].value) {
          return `$${(val / suffixes[i].value).toFixed(1)}${suffixes[i].symbol}`
        }
      }

      // fallback for small values
      return d3.format('$,')(val)
    }

    let infoObj = null
    let level = d.depth

    if (level === 1 && familyInfo[d.data.name]) {
      // family level
      const familyData = familyInfo[d.data.name]
      infoObj = {
        type: 'family',
        name: d.data.name,
        image: familyData.image,
        text: familyData.text,
        caption: familyData.caption,
        economicValue: d3.format('$.1s')(d.value).replace('G', 'B'),
        speciesCount: d.data.originalSpeciesCount || (d.children ? d.children.length : 0)
      }
    } else if (level === 2) {
      // species level
      const familyName = d.parent?.data?.name
      const silhouette = familyInfo[familyName]?.image

      let economicValue
      let lumpedSpeciesCount = 0
      if (d.data.name === 'Other') {
        const totalValue = d.data.children.reduce((sum, c) => sum + (c.value || 0), 0)
        economicValue = totalValue
        lumpedSpeciesCount = d.data.lumpedSpeciesCount || 0
      } else {
        economicValue = formatEconomicValue(d.value)
      }

      infoObj = {
        type: 'species',
        name: d.data.name,
        family: familyName,
        image: silhouette,
        economicValue,
        countryCount: d.children ? d.children.length : 0,
        lumpedSpeciesCount
      }
    }
    activeFamily.value = infoObj || defaultFamily
  }

  return svg.node()
}
</script>

<style lang="scss">
/* css for elements added/classed w/ d3 */
#circle-pack-svg {
  max-width: 100%;
  height: auto;
}
.circle-packing-text {
  font-size: 1.2rem;
  stroke: white;
  stroke-width: 2px;
  paint-order: stroke;
  stroke-linejoin: round;
  @media screen and (max-width: 600px) {
    font-size: 2rem;
  }
}
.fish-title {
  font-weight: 700;
}
.fish-value {
  font-weight: 400;
}
</style>
