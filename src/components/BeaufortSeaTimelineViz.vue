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
            <button id="button-0" @click="updateChart">0</button>
            <button id="button-100" @click="updateChart">100</button>
            <button id="button-1200" @click="updateChart">1200</button>
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
    const chartDecade = ref(null);
    const chartWidth = 1152;
    const simulation = ref(null)

    // Declare behavior on mounted
    // functions called here
    onMounted(async () => {

        // Load the data
        const data = await d3.csv(publicPath + 'beaufort_species_abundance.csv');

        chartData.value = data
        chartDecade.value = '500'

        initChart({
            width: chartWidth, // outer width, in pixels
        });

        // build chart
        drawChart(chartData.value, {
            width: chartWidth,
            height: chartWidth,
            decade: chartDecade.value
        })
        // BubbleChart(chartData.value, {
        //     label: d => [...d.id.split(".").pop().split(/(?=[A-Z][a-z])/g), d.value.toLocaleString("en")].join("\n"),
        //     value: d => d.value,
        //     id: d => d.id,
        //     group: d => d.id.split(".")[1],
        //     title: d => `${d.id}\n${d.value.toLocaleString("en")}`,
        //     link: d => `https://github.com/prefuse/Flare/blob/master/flare/src/${d.id.replace(/\./g, "/")}.as`,
        //     width: chartWidth,
        //     year: chartYear.value
        // })
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
    function drawChart(data, {
        width = 200,
        height = 200,
        decade = 200
    }) {
        
        data = data.filter(d => d.decade === decade);
        
        // Set radius based on data value
        const sizeScale = d3.scaleSqrt()
            .domain([
                d3.min(data, (d) => parseFloat(d.pct_abundance)),
                d3.max(data, (d) => parseFloat(d.pct_abundance))
            ])
            .range([4, 200]);

        const nodes = data.map((d) => ({
            ...d,
            radius: sizeScale(parseFloat(d.pct_abundance)),
            x: width / 2,
            y: height / 2
        }));
       
        // set up nodes
        let nodeGroups = svg.value.selectAll('.nodes')
            .data(nodes, d => d.species)
            // .enter()
            // .append("circle")
            // .attr("r", (d) => d.radius)
            // .attr("cx", width / 2)
            // .attr("cy", height / 2)
            // .attr("fill", (d) => d.hexcode)
            // .attr("stroke", "#000000")
            // .attr("stroke-width", 0.5);
        console.log(nodeGroups)

        const oldNodeGroups = nodeGroups.exit()

        // oldNodeGroups.selectAll("circle")
            // .attr("r", 0); //radius to 0

        // Remove old nodes
        oldNodeGroups.transition(getExitTransition()).remove()

        // Append new nodes
        const newNodeGroups = nodeGroups.enter().append("g")
            .attr("class", "nodes")
            .attr("id", d => d.id)

        newNodeGroups.append("circle")
            .attr("stroke", "#000000")
            .attr("stroke-width", 0.5)
            .attr("fill", d => d.hexcode)
            .attr("cx", width / 2)
            .attr("cy", height / 2)
            .attr("r", d => d.radius) //instantiate w/ radius = 0

        // const newNodeGroups = nodeGroups
        //     .enter()
        //     .append("circle")
        //     .attr("stroke", "#000000")
        //     .attr("stroke-width", 0.5)
        //     .attr("fill", d => d.hexcode)
        //     .attr("cx", width / 2)
        //     .attr("cy", height / 2)
        //     .attr("r", 0)

        //update nodeGroups to include new nodes
        nodeGroups = newNodeGroups.merge(nodeGroups)

        const nodeGroupCircle = nodeGroups.select("circle")
        console.log(nodeGroupCircle)

        nodeGroupCircle
            .transition(getUpdateTransition())
            .attr("r", d => d.radius)
            // .attr("cx", (d, i) => d.x + i*10)
            // .attr("cy", (d, i) => d.y + i*10);

        // nodeGroups
        //     .transition(getUpdateTransition())
        //     .attr("r", d => d.radius)

        // const node = svg.value
        //     .selectAll("circle")
        //     .data(nodes)
        //     .enter()
        //     .append("circle")
        //     .attr("r", (d) => d.radius)
        //     .attr("fill", (d) => d.color);

                function ticked() {
            nodeGroupCircle
            // nodeGroups
                .transition(getUpdateTransition())
                .attr("cx", (d) => d.x)
                .attr("cy", (d) => d.y);
        }


        simulation.value = d3.forceSimulation();
        simulation.value
            .nodes(nodes)
            .force("x", d3.forceX(width / 2).strength(0.05))
            .force("y", d3.forceY(height / 2).strength(0.05))
            .force(
                "collide",
                d3.forceCollide()
                    .radius((d) => d.radius + 2)
                    .iterations(1)
            )
            .force('charge', d3.forceManyBody().strength(0))
            .on("tick", ticked);



        simulation.value
            .on("tick", ticked);
        // if (simulation) {
        //     simulation.force("x", d3.forceX(width / 2).strength(0.1));
        //     simulation.alpha(1).restart();
        // }
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

        const r = d3.scaleSqrt(
            [100, 1000],
            [1, Math.sqrt()]
        )

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
        const simulation = d3.forceSimulation();
        simulation
            .nodes(root.leaves())
            .force("x", d3.forceX(width / 2).strength(0.1))
            .force("y", d3.forceY(height / 2).strength(0.1))
            .force(
                "collide",
                d3
                    .forceCollide()
                    .radius((d) => d.r)
                    .iterations(1)
            )
            .on("tick", ticked);

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

        function ticked() {
            leafGroupCircle
                .transition(getUpdateTransition())
                .attr("cx", (d) => d.x)
                .attr("cy", (d) => d.y);
        }

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

        // return Object.assign(svg.value.node(), {scales: {color}});
    }

    // define transitions
    function getUpdateTransition() {
      return d3.transition()
        .duration(100)
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

        chartDecade.value = e.target.id.split("-")[1]

        drawChart(chartData.value, {
            width: chartWidth,
            height: chartWidth,
            decade: chartDecade.value
        })

        // build chart
        // BubbleChart(chartData.value , {
        //     label: d => [...d.id.split(".").pop().split(/(?=[A-Z][a-z])/g), d.value.toLocaleString("en")].join("\n"),
        //     value: d => d.value,
        //     id: d => d.id,
        //     group: d => d.id.split(".")[1],
        //     title: d => `${d.id}\n${d.value.toLocaleString("en")}`,
        //     link: d => `https://github.com/prefuse/Flare/blob/master/flare/src/${d.id.replace(/\./g, "/")}.as`,
        //     width: 1152,
        //     year: chartYear.value
        // })
    }
</script>

<style>
</style>