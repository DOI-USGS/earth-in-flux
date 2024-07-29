<template>
    <!---VizSection-->
    <VizSection
        :figures="true"
        :fig-caption="false"
    >
        <!-- HEADING -->
        <!-- FIGURES -->
        <template #aboveExplanation>
        </template>
        <template #figures>
            <div class="single maxWidth">
                <figure>
                    <img src="https://labs.waterdata.usgs.gov/visualizations/images/BeaufortSeaTimeline.png">
                </figure>
            </div>
        </template>
        <!-- FIGURE CAPTION -->
        <template #figureCaption>
        </template>
        <!-- EXPLANATION -->
        <template #belowExplanation>
        </template>
    </VizSection>
    <!---VizSection-->
    <VizSection
        :figures="true"
        :fig-caption="false"
    >
        <!-- HEADING -->
        <template #heading>
        </template>
        <!-- FIGURES -->
        <template #aboveExplanation>
            <button id="button-2010" @click="updateChart">2010</button>
            <button id="button-2020" @click="updateChart">2020</button>
        </template>
        <template #figures>
            <div class="chart-container single" ref="chart">
            </div>
        </template>
        <!-- FIGURE CAPTION -->
        <template #figureCaption>
        </template>
        <!-- EXPLANATION -->
        <template #belowExplanation>
        </template>
    </VizSection>
</template>

<script setup>
    import { onMounted, ref } from "vue";
    import * as d3 from 'd3';
    import VizSection from '@/components/VizSection.vue';

    // define props
    defineProps({
        text: { type: Object }
    })

    // global variables
    const publicPath = import.meta.env.BASE_URL;
    const chart = ref(null);
    const svg = ref(null);
    const chartData = ref(null);
    const chartYear = ref(null);
    const chartWidth = 1152;

    // Declare behavior on mounted
    // functions called here
    onMounted(async () => {

        // Load the data
        const data = await d3.csv(publicPath + 'flare.csv');

        chartData.value = data.filter(d => d.value !== null) // just the leaves

        chartYear.value = '2010'

        

        initChart({
            width: chartWidth, // outer width, in pixels
        });

        // build chart
        BubbleChart(chartData.value, {
            label: d => [...d.id.split(".").pop().split(/(?=[A-Z][a-z])/g), d.value.toLocaleString("en")].join("\n"),
            value: d => d.value,
            id: d => d.id,
            group: d => d.id.split(".")[1],
            title: d => `${d.id}\n${d.value.toLocaleString("en")}`,
            link: d => `https://github.com/prefuse/Flare/blob/master/flare/src/${d.id.replace(/\./g, "/")}.as`,
            width: chartWidth,
            year: chartYear.value
        })
    });

    function initChart({
        width = 640, // outer width, in pixels
        height = width, // outer height, in pixels
        margin = 1, // default margins
        marginTop = margin, // top margin, in pixels
        marginLeft = margin // left margin, in pixels
    }) {
        svg.value = d3.select(chart.value)
            .append("svg")
            .attr("id", "circle-pack-svg")
            .attr("width", width)
            .attr("height", height)
            .attr("viewBox", [-marginLeft, -marginTop, width, height])
            .attr("style", "max-width: 100%; height: auto; height: intrinsic;")
            .attr("fill", "currentColor")
            .attr("font-size", 10)
            .attr("font-family", "sans-serif")
            .attr("text-anchor", "middle");
    }

    // Copyright 2021-2023 Observable, Inc.
    // Released under the ISC license.
    // https://observablehq.com/@d3/bubble-chart
    function BubbleChart(data, {
        name = ([x]) => x, // alias for label
        label = name, // given d in data, returns text to display on the bubble
        value = ([, y]) => y, // given d in data, returns a quantitative size
        id = '',
        group, // given d in data, returns a categorical value for color
        title, // given d in data, returns text to show on hover
        link, // given a node d, its link (if any)
        linkTarget = "_blank", // the target attribute for links, if any
        width = 640, // outer width, in pixels
        height = width, // outer height, in pixels
        padding = 3, // padding between circles
        margin = 1, // default margins
        marginTop = margin, // top margin, in pixels
        marginRight = margin, // right margin, in pixels
        marginBottom = margin, // bottom margin, in pixels
        marginLeft = margin, // left margin, in pixels
        groups, // array of group names (the domain of the color scale)
        colors = d3.schemeTableau10, // an array of colors (for groups)
        fill = "#ccc", // a static fill color, if no group channel is specified
        fillOpacity = 0.7, // the fill opacity of the bubbles
        stroke, // a static stroke around the bubbles
        strokeWidth, // the stroke width around the bubbles, if any
        strokeOpacity, // the stroke opacity around the bubbles, if any
        year // the year for which to draw data
        } = {}) {
        // Compute the values.
        data = data.filter(d => d.year === year);
        const D = d3.map(data, d => d);
        const V = d3.map(data, value);
        const G = group == null ? null : d3.map(data, group);
        const I = d3.range(V.length).filter(i => V[i] > 0);
        const ID = d3.map(data, id)

        // Unique the groups.
        if (G && groups === undefined) groups = I.map(i => G[i]);
        groups = G && new d3.InternSet(groups);

        // Construct scales.
        const color = G && d3.scaleOrdinal(groups, colors);

        // Compute labels and titles.
        const L = label == null ? null : d3.map(data, label);
        const T = title === undefined ? L : title == null ? null : d3.map(data, title);

        // Compute layout: create a 1-deep hierarchy, and pack it.
        const root = d3.pack()
            .size([width - marginLeft - marginRight, height - marginTop - marginBottom])
            .padding(padding)
            (d3.hierarchy({children: I})
            .sum(i => V[i]));

        // const svg = d3.select(chart.value)
        //     .append("svg")
        //     .attr("id", "circle-pack-svg")
        //     .attr("width", width)
        //     .attr("height", height)
        //     .attr("viewBox", [-marginLeft, -marginTop, width, height])
        //     .attr("style", "max-width: 100%; height: auto; height: intrinsic;")
        //     .attr("fill", "currentColor")
        //     .attr("font-size", 10)
        //     .attr("font-family", "sans-serif")
        //     .attr("text-anchor", "middle");

        // const leaf = svg.selectAll("g")
        //     .data(root.leaves(), d => ID[d.data])
        //     .join("g")
        //     .attr("transform", d => `translate(${d.x},${d.y})`);

        let leafGroups = svg.value.selectAll('.leaves')
            .data(root.leaves(), d => ID[d.data])

        const oldLeafGroups = leafGroups.exit()

        oldLeafGroups.selectAll("circle")
            .attr("r", 0); //radius to 0

        // Remove old leaves
        oldLeafGroups.transition(getExitTransition()).remove()

        const newLeafGroups = leafGroups.enter().append("g")
            .attr("class", "leaves")
            .attr("id", d => ID[d.data])
            .attr("transform", d => `translate(${d.x},${d.y})`);

        newLeafGroups.append("circle")
            .attr("id", d => ID[d.data])
            .attr("stroke", stroke)
            .attr("stroke-width", strokeWidth)
            .attr("stroke-opacity", strokeOpacity)
            .attr("fill", G ? d => color(G[d.data]) : fill == null ? "none" : fill)
            .attr("fill-opacity", fillOpacity)
            .attr("r", 0); //instantiate w/ radius = 0

        // update rectGroups to include new points
        leafGroups = newLeafGroups.merge(leafGroups)

        leafGroups
            .transition(getUpdateTransition())
            .attr("transform", d => `translate(${d.x},${d.y})`);

        const leafGroupCircle = leafGroups.select("circle")

        leafGroupCircle
            .transition(getUpdateTransition())
            .attr("r", d => d.r);

        if (T) leafGroups.append("title")
            .text(d => T[d.data]);

        // if (L) {
        //     // A unique identifier for clip paths (to avoid conflicts).
        //     const uid = `O-${Math.random().toString(16).slice(2)}`;

        // leafGroups.append("clipPath")
        //         .attr("id", d => `${uid}-clip-${d.data}`)
        //     .append("circle")
        //         .attr("r", d => d.r);

        // leafGroups.append("text")
        //         .attr("clip-path", d => `url(${new URL(`#${uid}-clip-${d.data}`, location)})`)
        //     .selectAll("tspan")
        //     .data(d => `${L[d.data]}`.split(/\n/g))
        //     .join("tspan")
        //         .attr("x", 0)
        //         .attr("y", (d, i, D) => `${i - D.length / 2 + 0.85}em`)
        //         .attr("fill-opacity", (d, i, D) => i === D.length - 1 ? 0.7 : null)
        //         .text(d => d);
        // }

        // leaf.append("circle")
        //     .attr("id", d => ID[d.data])
        //     .attr("stroke", stroke)
        //     .attr("stroke-width", strokeWidth)
        //     .attr("stroke-opacity", strokeOpacity)
        //     .attr("fill", G ? d => color(G[d.data]) : fill == null ? "none" : fill)
        //     .attr("fill-opacity", fillOpacity)
        //     .attr("r", d => d.r);

        // if (T) leaf.append("title")
        //     .text(d => T[d.data]);

        // if (L) {
        //     // A unique identifier for clip paths (to avoid conflicts).
        //     const uid = `O-${Math.random().toString(16).slice(2)}`;

        //     leaf.append("clipPath")
        //         .attr("id", d => `${uid}-clip-${d.data}`)
        //     .append("circle")
        //         .attr("r", d => d.r);

        //     leaf.append("text")
        //         .attr("clip-path", d => `url(${new URL(`#${uid}-clip-${d.data}`, location)})`)
        //     .selectAll("tspan")
        //     .data(d => `${L[d.data]}`.split(/\n/g))
        //     .join("tspan")
        //         .attr("x", 0)
        //         .attr("y", (d, i, D) => `${i - D.length / 2 + 0.85}em`)
        //         .attr("fill-opacity", (d, i, D) => i === D.length - 1 ? 0.7 : null)
        //         .text(d => d);
        // }

        return Object.assign(svg.value.node(), {scales: {color}});
    }

    // define transitions
    function getUpdateTransition() {
      return d3.transition()
        .duration(500)
        .ease(d3.easeCubicInOut)
    }
    function getExitTransition() {
      return d3.transition()
        .duration(500)
        .ease(d3.easeCubicInOut)
    }
    function updateChart(e) {
        console.log('clicked!')
        console.log(e.target.id.split("-")[1])

        chartYear.value = e.target.id.split("-")[1]

        // build chart
        BubbleChart(chartData.value , {
            label: d => [...d.id.split(".").pop().split(/(?=[A-Z][a-z])/g), d.value.toLocaleString("en")].join("\n"),
            value: d => d.value,
            id: d => d.id,
            group: d => d.id.split(".")[1],
            title: d => `${d.id}\n${d.value.toLocaleString("en")}`,
            link: d => `https://github.com/prefuse/Flare/blob/master/flare/src/${d.id.replace(/\./g, "/")}.as`,
            width: 1152,
            year: chartYear.value
        })
    }
</script>

<style>
</style>