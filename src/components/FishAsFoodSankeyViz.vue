<template>
  <!---VizSection-->
  <VizSection :figures="true" :fig-caption="false">
    <!-- HEADING -->
    <template #heading> </template>
    <!-- FIGURES -->
    <template #aboveExplanation>
      <p v-html="text.paragraph1" />
      <p class="annotation">{{ currentText }}</p>
    </template>
    <template #figures>
      <div class="chart-container single" ref="chart"></div>
    </template>
    <!-- FIGURE CAPTION -->
    <template #figureCaption> </template>
    <!-- EXPLANATION -->
    <template #belowExplanation> </template>
  </VizSection>
</template>

<script setup>
import { onMounted, ref } from 'vue'
import * as d3 from 'd3'
import * as d3Sankey from 'd3-sankey'
import VizSection from '@/components/VizSection.vue'
import { isMobile } from 'mobile-device-detect'

defineProps({
  text: { type: Object }
})

const publicPath = import.meta.env.BASE_URL
const chart = ref(null)
const nodeAlign = 'justify'
const currentText = ref('')
const mobileView = isMobile

// Family color palette
const colors = ['#2b2e3c', '#5b7083', '#cc5b4d', '#d09a47', '#628c8c']

onMounted(async () => {
  const fontSize = mobileView ? 11 : 13
  const chartWidth = mobileView ? chart.value.clientWidth : 800
  const chartHeight = mobileView ? 1600 : 1200

  const rawLinks = await d3.csv(publicPath + 'fish_as_food_harvest.csv')

  // Get top-level families
  const allSources = new Set(rawLinks.map((d) => d.source))
  const allTargets = new Set(rawLinks.map((d) => d.target))
  const families = [...allSources].filter((src) => !allTargets.has(src))

  const colorScale = d3.scaleOrdinal().domain(families).range(colors)

  // Trace every node back to their family
  const parentOf = new Map()
  rawLinks.forEach((d) => parentOf.set(d.target, d.source))

  const familyGroupMap = new Map()

  function findFamily(node) {
    if (families.includes(node)) {
      familyGroupMap.set(node, node)
      return node
    }
    const parent = parentOf.get(node)
    if (!parent) return null
    const fam = findFamily(parent)
    familyGroupMap.set(node, fam)
    return fam
  }

  rawLinks.forEach((d) => {
    findFamily(d.source)
    findFamily(d.target)
  })

  // Find rhe dominant fish family harvested in each country and then color the nodes accordingly
  const countryFamilyTotals = new Map() // create map to track total harvest vals

  // loop through each link in data
  rawLinks.forEach((d) => {
    // check if the target is a country by confirming it's *not* in the list of sources
    const isCountry = !allSources.has(d.target)

    if (isCountry) {
      const country = d.target // destination country
      const species = d.source // species being harvested
      const family = findFamily(species) // determine the family this species belongs to
      const value = +d.value // convert  harvest value to number

      // if this country hasn't been seen yet, initialize its entry in the map with an empty object
      if (!countryFamilyTotals.has(country)) {
        countryFamilyTotals.set(country, {})
      }

      //  get the totals object for this country
      const totals = countryFamilyTotals.get(country)
      // if family hasn't been counted for this country, initialize it to 0
      if (!totals[family]) {
        totals[family] = 0
      }
      // Add  current value to the running total for this family in this country
      totals[family] += value
    }
  })

  // set color for each country based on dominant family
  const countryColorMap = new Map()

  countryFamilyTotals.forEach((famObj, country) => {
    const [topFamily] = Object.entries(famObj).sort((a, b) => b[1] - a[1])[0]
    countryColorMap.set(country, colorScale(topFamily))
  })

  // Pass everything into SankeyChart
  SankeyChart(
    { links: rawLinks },
    {
      nodeGroup: (d) => familyGroupMap.get(d.id),
      nodeGroups: families,
      nodeAlign,
      format: d3.format(',.0f'),
      width: chartWidth,
      height: chartHeight,
      fontSize,
      colors,
      nodeStroke: 'none',
      nodeLabelPadding: mobileView ? 2 : 6,
      linkColor: (d) => colorScale(familyGroupMap.get(d.source.id)),
      countryColorMap // pass to chart
    }
  )
})

// https://observablehq.com/@d3/sankey-component
// Copyright 2021-2023 Observable, Inc.
// Released under the ISC license.
// https://observablehq.com/@d3/sankey-diagram
function SankeyChart(
  {
    nodes, // an iterable of node objects (typically [{id}, …]); implied by links if missing
    links // an iterable of link objects (typically [{source, target}, …])
  },
  {
    format = ',', // a function or format specifier for values in titles
    align = 'justify', // convenience shorthand for nodeAlign
    nodeId = (d) => d.id, // given d in nodes, returns a unique identifier (string)
    nodeGroup, // given d in nodes, returns an (ordinal) value for color
    nodeGroups, // an array of ordinal values representing the node groups
    nodeLabel, // given d in (computed) nodes, text to label the associated rect
    nodeTitle = (d) => `${d.id}\n${format(d.value)} kg`, // given d in (computed) nodes, hover text
    nodeAlign = align, // Sankey node alignment strategy: left, right, justify, center
    nodeSort, // comparator function to order nodes
    nodeWidth = 15, // width of node rects
    nodePadding = 10, // vertical separation between adjacent nodes
    nodeLabelPadding = 6, // horizontal separation between node and label
    nodeStroke = 'none', // stroke around node rects
    nodeStrokeWidth, // width of stroke around node rects, in pixels
    nodeStrokeOpacity, // opacity of stroke around node rects
    nodeStrokeLinejoin, // line join for stroke around node rects
    linkSource = ({ source }) => source, // given d in links, returns a node identifier string
    linkTarget = ({ target }) => target, // given d in links, returns a node identifier string
    linkValue = ({ value }) => value, // given d in links, returns the quantitative value
    linkPath = d3Sankey.sankeyLinkHorizontal(), // given d in (computed) links, returns the SVG path
    linkTitle = (d) => `${d.source.id} → ${d.target.id}\n${format(d.value)} kg`, // given d in (computed) links
    linkColor = 'source-target', // source, target, source-target, or static color
    linkStrokeOpacity = 0.5, // link stroke opacity
    linkMixBlendMode = 'multiply', // link blending mode
    width = 640, // outer width, in pixels
    height = 400, // outer height, in pixels
    marginTop = 5, // top margin, in pixels
    marginRight = 1, // right margin, in pixels
    marginBottom = 5, // bottom margin, in pixels
    marginLeft = 1, // left margin, in pixels
    countryColorMap,
    fontSize = 12
  } = {}
) {
  // Convert nodeAlign from a name to a function (since d3-sankey is not part of core d3).
  if (typeof nodeAlign !== 'function')
    nodeAlign =
      {
        left: d3Sankey.sankeyLeft,
        right: d3Sankey.sankeyRight,
        center: d3Sankey.sankeyCenter
      }[nodeAlign] ?? d3Sankey.sankeyJustify

  // Compute values.
  const LS = d3.map(links, linkSource).map(intern)
  const LT = d3.map(links, linkTarget).map(intern)
  const LV = d3.map(links, linkValue)
  if (nodes === undefined) nodes = Array.from(d3.union(LS, LT), (id) => ({ id }))
  const N = d3.map(nodes, nodeId).map(intern)
  const G = nodeGroup == null ? null : d3.map(nodes, nodeGroup).map(intern)

  // Replace the input nodes and links with mutable objects for the simulation.
  nodes = d3.map(nodes, (_, i) => ({ id: N[i] }))
  links = d3.map(links, (_, i) => ({ source: LS[i], target: LT[i], value: LV[i] }))

  // Ignore a group-based linkColor option if no groups are specified.
  if (!G && ['source', 'target', 'source-target'].includes(linkColor)) linkColor = 'currentColor'

  // Compute default domains.
  if (G && nodeGroups === undefined) nodeGroups = G

  // Construct the scales.
  const color = nodeGroup == null ? null : d3.scaleOrdinal(nodeGroups, colors)

  // Compute the Sankey layout.
  d3Sankey
    .sankey()
    .nodeId(({ index: i }) => N[i])
    .nodeAlign(nodeAlign)
    .nodeWidth(nodeWidth)
    .nodePadding(nodePadding)
    .nodeSort(nodeSort)
    .extent([
      [marginLeft, marginTop],
      [width - marginRight, height - marginBottom]
    ])({ nodes, links })

  // Compute titles and labels using layout nodes, so as to access aggregate values.
  if (typeof format !== 'function') format = d3.format(format)
  const Tl = nodeLabel === undefined ? N : nodeLabel == null ? null : d3.map(nodes, nodeLabel)
  const Tt = nodeTitle == null ? null : d3.map(nodes, nodeTitle)
  const Lt = linkTitle == null ? null : d3.map(links, linkTitle)

  // A unique identifier for clip paths (to avoid conflicts).
  const uid = `O-${Math.random().toString(16).slice(2)}`

  const svg = d3
    .select(chart.value)
    .append('svg')
    .attr('id', 'sankey-svg')
    .attr('width', width)
    .attr('height', height)
    .attr('viewBox', [0, 0, width, height])
    .attr('style', 'max-width: 100%; height: auto; height: intrinsic;')

  const node = svg
    .append('g')
    .attr('stroke', nodeStroke)
    .attr('stroke-width', nodeStrokeWidth)
    .attr('stroke-opacity', nodeStrokeOpacity)
    .attr('stroke-linejoin', nodeStrokeLinejoin)
    .selectAll('rect')
    .data(nodes)
    .join('rect')
    .attr('x', (d) => d.x0)
    .attr('y', (d) => d.y0)
    .attr('height', (d) => d.y1 - d.y0)
    .attr('width', (d) => d.x1 - d.x0)
    .attr('fill', (d) => {
      if (typeof countryColorMap !== 'undefined' && d.depth === 2 && countryColorMap.has(d.id)) {
        return countryColorMap.get(d.id)
      }
      return color(G[d.index])
    })
    .on('mouseover', (event, d) => {
      currentText.value = Tt[d.index]
    })
    .on('mouseout', () => {
      currentText.value = ''
    })

  if (Tt) node.append('title').text(({ index: i }) => Tt[i])

  const link = svg
    .append('g')
    .attr('fill', 'none')
    .attr('stroke-opacity', linkStrokeOpacity)
    .selectAll('g')
    .data(links)
    .join('g')
    .style('mix-blend-mode', linkMixBlendMode)

  if (linkColor === 'source-target')
    link
      .append('linearGradient')
      .attr('id', (d) => `${uid}-link-${d.index}`)
      .attr('gradientUnits', 'userSpaceOnUse')
      .attr('x1', (d) => d.source.x1)
      .attr('x2', (d) => d.target.x0)
      .call((gradient) =>
        gradient
          .append('stop')
          .attr('offset', '0%')
          .attr('stop-color', ({ source: { index: i } }) => color(G[i]))
      )
      .call((gradient) =>
        gradient
          .append('stop')
          .attr('offset', '100%')
          .attr('stop-color', ({ target: { index: i } }) => color(G[i]))
      )

  link
    .append('path')
    .attr('d', linkPath)
    .attr(
      'stroke',
      typeof linkColor === 'function'
        ? (d) => linkColor(d)
        : linkColor === 'source-target'
          ? ({ index: i }) => `url(#${uid}-link-${i})`
          : linkColor === 'source'
            ? ({ source: { index: i } }) => color(G[i])
            : linkColor === 'target'
              ? ({ target: { index: i } }) => color(G[i])
              : linkColor
    )
    .attr('stroke-width', ({ width }) => Math.max(1, width))
    .call(Lt ? (path) => path.append('title').text(({ index: i }) => Lt[i]) : () => {})
    .on('mouseover', (event, d) => {
      currentText.value = Lt[d.index]
    })
    .on('mouseout', () => {
      currentText.value = ''
    })

  if (Tl)
    svg
      .append('g')
      .selectAll('text')
      .data(nodes)
      .join('text')
      .attr('x', (d) => (d.x0 < width / 2 ? d.x1 + nodeLabelPadding : d.x0 - nodeLabelPadding))
      .attr('y', (d) => (d.y1 + d.y0) / 2)
      .attr('dy', '0.35em')
      .attr('text-anchor', (d) => (d.x0 < width / 2 ? 'start' : 'end'))
      .attr('font-size', fontSize)
      .attr('fill', 'black')
      .attr('stroke', 'white')
      .attr('stroke-width', 2)
      .attr('paint-order', 'stroke')
      .attr('stroke-linejoin', 'round')
      .style('user-select', 'none')
      .text(({ index: i }) => Tl[i])

  function intern(value) {
    return value !== null && typeof value === 'object' ? value.valueOf() : value
  }

  return Object.assign(svg.node(), { scales: { color } })
}
</script>

<style scoped lang="scss">
.chart-container {
  margin-top: 4rem;

  @media (max-width: 700px) {
    overflow-x: auto;
    padding-bottom: 1rem;
  }

  svg {
    @media (max-width: 700px) {
      width: 100% !important;
      height: auto !important;
    }
  }
}
.annotation {
  margin-top: 3rem;
  font-style: italic;
  font-weight: 300;
  height: 1rem;
}
</style>
