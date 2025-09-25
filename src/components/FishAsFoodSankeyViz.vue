<template>
  <!---VizSection-->
  <VizSection :figures="true" :fig-caption="false">
    <!-- HEADING -->
    <template #heading> </template>
    <!-- FIGURES -->
    <template #aboveExplanation>
      <p v-html="text.paragraph1" />
      <div class="toggle-group">
        <ToggleSwitch 
            v-for="toggle, index in toggles"
            :key="index"
            v-model="toggle.value" 
            :label="toggle.label"
            :rightColor="toggle.color"
        />
      </div>
      <p class="annotation">{{ currentText }}</p>
    </template>
    <template #figures>
      <div class="chart-container single" ref="chart"></div>
    </template>
    <!-- EXPLANATION -->
    <template #belowExplanation>
      <p v-html="text.paragraph2" />
    </template>
  </VizSection>
</template>

<script setup>
import { onMounted, reactive, ref, watch } from 'vue'
import * as d3 from 'd3'
import * as d3Sankey from 'd3-sankey'
import VizSection from '@/components/VizSection.vue'
import ToggleSwitch from '@/components/ToggleSwitch.vue'
import { isMobile } from 'mobile-device-detect'

const props = defineProps({
  text: { type: Object }
})

// global variables
const publicPath = import.meta.env.BASE_URL
const chart = ref(null)
const nodeAlign = 'justify'
const currentText = ref('')
const mobileView = isMobile

// set up reactive layers object
const toggles = reactive(props.text.toggleData);

// Family color palette
const colors = {
  'Bagridae' : '#2b2e3c',
  'Centrarchidae' : '#5b7083',
  'Cyprinidae' : '#cc5b4d',
  'Percidae' : '#d09a47',
  'Salmonidae' : '#628c8c'
}

let rawLinks = null;
let filteredLinks = null;

onMounted(async () => {
  rawLinks = await d3.csv(publicPath + 'fish_as_food_harvest.csv', (d) => ({
    source: d.source,
    target: d.target,
    value: +d.value
  }))

  // Build filtered set of links that don't include data for China

  // Find species/families that ONLY link to China
  const allTargetsFromSource = new Map()
  rawLinks.forEach((d) => {
    if (!allTargetsFromSource.has(d.source)) {
      allTargetsFromSource.set(d.source, new Set())
    }
    allTargetsFromSource.get(d.source).add(d.target)
  })
  const chinaOnlySources = new Set()
  allTargetsFromSource.forEach((targets, source) => {
    if (targets.size === 1 && targets.has('China')) {
      chinaOnlySources.add(source)
    }
  })

  // Build Map of species + values associated with China
  const chinaLinks = rawLinks.filter((d) => d.target == 'China')
  const chinaSourceValues = new Map()
  chinaLinks.forEach((d) => {
    chinaSourceValues.set(d.source, d.value)
  })

  // Generate filtered links
  filteredLinks = rawLinks
    // First filter out links with China as a target
    .filter((d) => d.target !== 'China')
    // Filter out links with sources that are China-only sources
    .filter((d) => !chinaOnlySources.has(d.source))
    // Filter out links with targets that are China-only sources
    .filter((d) => !chinaOnlySources.has(d.target))
    // make a deep copy, so that we can modify this object without affecting rawLinks
    .map(d => JSON.parse(JSON.stringify(d)))
    // For source links with targets that are China sources and sources for other countries
    // update values to remove value totals associated with China
    .map((d) => {
      if (chinaSourceValues.has(d.target)) {
        return {...d, value: d.value-= chinaSourceValues.get(d.target)}
      } else {
        return d
      }
    })

  // draw sankey 
  const showChina = toggles.showChina.value;
  const sortByHarvest = toggles.sortHarvest.value;
  const dataset = showChina ? rawLinks : filteredLinks
  drawSankey(dataset, showChina, sortByHarvest)
})

watch(
  [toggles],
  ([togglesNewVal]) => {
    const showChina = togglesNewVal.showChina.value;
    const sortByHarvest = togglesNewVal.sortHarvest.value;
    const dataset = showChina ? rawLinks : filteredLinks
    drawSankey(dataset, showChina, sortByHarvest)
  }
)

// Tie colors to families
function makeColorScale(categories, colors) {
  return d3.scaleOrdinal()
    .domain(categories)
    .range(categories.map(category => colors[category]));
}

function drawSankey(links, showChina, sortByHarvest) {
  d3.select(chart.value).selectAll('*').remove()

  if (!links || links.length === 0) return

  const validLinks = links.filter((d) => d.value && !isNaN(d.value))

  const allSources = new Set(validLinks.map((d) => d.source))
  const allTargets = new Set(validLinks.map((d) => d.target))
  const families = [...allSources].filter((src) => !allTargets.has(src))

  const colorScale = makeColorScale(families, colors)

  const parentOf = new Map()
  validLinks.forEach((d) => parentOf.set(d.target, d.source))

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

  validLinks.forEach((d) => {
    findFamily(d.source)
    findFamily(d.target)
  })

  const nodes = Array.from(
    new Set([...validLinks.map((d) => d.source), ...validLinks.map((d) => d.target)])
  ).map((id) => ({ id }))

  const activeNodes = nodes.filter((node) =>
    validLinks.some((link) => link.source === node.id || link.target === node.id)
  )

  // Country colors
  const countryFamilyTotals = new Map()
  validLinks.forEach((d) => {
    const isCountry = !allSources.has(d.target)
    if (isCountry) {
      const country = d.target
      const species = d.source
      const family = findFamily(species)
      const value = +d.value
      if (!countryFamilyTotals.has(country)) countryFamilyTotals.set(country, {})
      const totals = countryFamilyTotals.get(country)
      totals[family] = (totals[family] || 0) + value
    }
  })

  const countryColorMap = new Map()
  countryFamilyTotals.forEach((famObj, country) => {
    const [topFamily] = Object.entries(famObj).sort((a, b) => b[1] - a[1])[0]
    countryColorMap.set(country, colorScale(topFamily))
  })

  const desktopHeight = showChina ? 1200 : 900;
  const mobileHeight = showChina ? 1600 : 1200;
  SankeyChart(
    { nodes: activeNodes, links: validLinks },
    {
      nodeGroup: (d) => familyGroupMap.get(d.id),
      nodeGroups: families,
      nodeAlign,
      nodeSort: sortByHarvest ? (a,b) => d3.descending(a.value, b.value) : undefined,
      format: d3.format(',.0f'),
      width: mobileView ? chart.value.clientWidth : 800,
      height: mobileView ? mobileHeight : desktopHeight,
      fontSize: mobileView ? 11 : 13,
      colors,
      nodeStroke: 'none',
      nodeLabelPadding: mobileView ? 2 : 6,
      linkColor: (d) => colorScale(familyGroupMap.get(d.source.id)),
      countryColorMap
    }
  )
}

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
  const color = nodeGroup == null ? null : makeColorScale(nodeGroups, colors)

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
  margin-bottom: 4rem;
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
