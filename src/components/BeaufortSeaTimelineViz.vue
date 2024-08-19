<template>
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
            <button v-for="decade in chartDecades" :key="decade" :id="`button-${decade}`" @click="updateChart" :class="chartDecade === decade ? 'active' : 'inactive'">
                {{ decade }}
            </button>
        </template>
        <template #figures>
            <div id="chart-container"></div>
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
    const dataFile = 'beaufort_species_abundance.csv';
    const data = ref(null);
    let chartDimensions;
    const chartTitle = 'Title of chart';
    let chartBounds;
    let chartDecades = ref(null);
    const chartDecade = ref(null);
    const simulation = ref(null)

    // Declare behavior on mounted
    // functions called here
    onMounted(async () => {
        try {
            await loadDatasets(dataFile);
            if (data.value.length > 0) {
                // pull array of unique decades
                chartDecades.value = Array.from(new Set(data.value.map(d => d.decade)))
                // set initial decade view
                chartDecade.value = chartDecades.value[0]
                
                // initialize svg and chart bounds
                initBubbleChart({
                    width: 900, // outer width, in pixels
                });

                // build chart
                drawBubbleChart(data.value, {
                    decade: chartDecade.value
                })
            } else {
                console.error('Error loading data');
            }
        } catch (error) {
            console.error('Error during component mounting', error);
        }
    });

    async function loadDatasets(file) {
        try {
            data.value = await loadData(file);
            console.log('data in');
        } catch (error) {
            console.error('Error loading datasets', error);
        }
    }

    async function loadData(fileName) {
        try {
            const data = await d3.csv(publicPath + fileName);
            return data;
        } catch (error) {
            console.error(`Error loading data from ${fileName}`, error);
            return [];
        }
    }

    function initBubbleChart({
        width = 640, // outer width, in pixels
        height = width, // outer height, in pixels
        margin = 1, // default margins
        marginTop = margin, // top margin, in pixels
        marginBottom = margin, // left margin, in pixels
        marginLeft = margin, // left margin, in pixels
        marginRight = margin // right margin, in pixels
    }) {
        // set up global chart dimensions
        chartDimensions = {
            width,
            height,
            margin: {
                top: marginTop,
                right: marginRight,
                bottom: marginBottom,
                left: marginLeft
            }
        }
        chartDimensions.boundedWidth = chartDimensions.width - chartDimensions.margin.left - chartDimensions.margin.right
        chartDimensions.boundedHeight = chartDimensions.height - chartDimensions.margin.top - chartDimensions.margin.bottom

        // draw canvas for chart
        const chartSVG = d3.select("#chart-container")
            .append("svg")
                .attr("viewBox", [0, 0, (chartDimensions.width), (chartDimensions.height)].join(' '))
                .attr("width", "100%")
                .attr("height", "100%")
                .attr("id", "chart-svg")

        // assign role for accessibility
        chartSVG.attr("role", "figure")
            .append("title")
            .text(chartTitle)

        // Add group for bounds
        chartBounds = chartSVG.append("g")
            .attr("id", "chart-bounds")
            .style("transform", `translate(${
                chartDimensions.margin.left
            }px, ${
                chartDimensions.margin.top
            }px)`)

        chartBounds.append("g")
            .attr("class", "nodes")
    }
    function drawBubbleChart(data, {
        decade = 200
    }) {
        // Set radius based on data values across all decades
        const sizeScale = d3.scaleSqrt()
            .domain([
                d3.min(data, (d) => parseFloat(d.pct_abundance)),
                d3.max(data, (d) => parseFloat(d.pct_abundance))
            ])
            .range([4, 200]);

        // filter data to current decade
        data = data.filter(d => d.decade === decade);
        
        const nodes = data.map((d) => ({
            ...d,
            radius: sizeScale(parseFloat(d.pct_abundance)),
            x: chartDimensions.boundedWidth / 2,
            y: chartDimensions.boundedHeight / 2
        }));
       
        // set up nodes
        let nodeGroups = chartBounds.selectAll('.nodes')
            .selectAll(".node")
            .data(nodes, d => d.species_id)

        const oldNodeGroups = nodeGroups.exit()

        oldNodeGroups.selectAll("circle")
            .attr("r", 0); //radius to 0

        // Remove old nodes
        oldNodeGroups.transition(getExitTransition()).remove()

        // Append new nodes
        const newNodeGroups = nodeGroups.enter().append("g")
            .attr("class", "node")
            .attr("id", d => "group_" + d.species_id)

        newNodeGroups.append("circle")
            .attr("id", d => d.species_id)
            .attr("stroke", "#000000")
            .attr("stroke-width", 0.5)
            .attr("fill", d => d.hexcode)
            .attr("r", 0) //instantiate w/ radius = 0

        //update nodeGroups to include new nodes
        nodeGroups = newNodeGroups.merge(nodeGroups)

        const nodeGroupCircle = nodeGroups.select("circle")
        
        nodeGroupCircle
            .transition(getUpdateTransition())
            .attr("r", d => d.radius)

        function ticked() {
            nodeGroupCircle
                // .transition(getUpdateTransition()) // BREAKS d3 force
                // .attr("r", d => d.radius)
                .attr("cx", (d) => d.x)
                .attr("cy", (d) => d.y);
        }

        // set up d3 force simulation
        if (simulation.value) {
            // simulation.value.stop()
            simulation.value
                .nodes(nodes)
                .alpha(0.9)
                .restart()
                .force("x", d3.forceX(chartDimensions.boundedWidth / 2).strength(0.05))
                .force("y", d3.forceY(chartDimensions.boundedHeight / 2).strength(0.05))
                .force("center", d3.forceCenter(chartDimensions.boundedWidth / 2, chartDimensions.boundedHeight / 2))
                .force(
                    "collide",
                    d3.forceCollide()
                        .radius((d) => d.radius + 2)
                        .iterations(1)
                )
                .force('charge', d3.forceManyBody().strength(0))
                // .alphaMin(0.01)
                // .alpha(0)
                // .velocityDecay(0.9)
                .on("tick", ticked);
        } else {
            simulation.value = d3.forceSimulation();
            simulation.value
                .nodes(nodes)
                .force("x", d3.forceX(chartDimensions.boundedWidth / 2).strength(0.05))
                .force("y", d3.forceY(chartDimensions.boundedHeight / 2).strength(0.05))
                .force("center", d3.forceCenter(chartDimensions.boundedWidth / 2, chartDimensions.boundedHeight / 2))
                .force(
                    "collide",
                    d3.forceCollide()
                        .radius((d) => d.radius + 2)
                        .iterations(1)
                )
                .force('charge', d3.forceManyBody().strength(0))
                // .alphaMin(0.01)
                // .alphaDecay(0)
                // .velocityDecay(0.9)
                .on("tick", ticked);
        }
    }

    // define transitions
    function getUpdateTransition() {
      return d3.transition()
        .duration(2000)
        .ease(d3.easeCubicInOut)
    }
    function getExitTransition() {
      return d3.transition()
        .duration(500)
        .ease(d3.easeCubicInOut)
    }
    function updateChart(e) {

        const clickedID = e.target.id.split("-")[1]

        if (clickedID != chartDecade.value) {

            chartDecade.value = clickedID

            drawBubbleChart(data.value, {
                decade: chartDecade.value
            })

        }

        // build chart
        // BubbleChart(data.value , {
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
    .active {
        background-color: grey;
        color: white;
    }
</style>