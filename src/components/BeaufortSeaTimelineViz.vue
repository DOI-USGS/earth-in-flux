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
            <button v-for="decade in chartDecades" :key="decade" :id="`button-${decade}`" @click="updateChart" :class="chartDecade === decade ? 'active' : 'inactive'">
                {{ decade }}
            </button>
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
    const chartWidth = 900;
    const simulation = ref(null)


    const chartDecades = ['0', '500', '1000', '1500', '2000']

    // Declare behavior on mounted
    // functions called here
    onMounted(async () => {

        // Load the data
        const data = await d3.csv(publicPath + 'beaufort_species_abundance.csv');

        chartData.value = data
        chartDecade.value = chartDecades[0]

        initChart({
            width: chartWidth, // outer width, in pixels
        });

        // build chart
        drawChart(chartData.value, {
            width: chartWidth,
            height: chartWidth,
            decade: chartDecade.value
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

        svg.value.append("g")
            .attr("class", "nodes")
    }
    function drawChart(data, {
        width = 200,
        height = 200,
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
            x: width / 2,
            y: height / 2
        }));
       
        // set up nodes
        let nodeGroups = svg.value.selectAll('.nodes')
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
                .force("x", d3.forceX(width / 2).strength(0.05))
                .force("y", d3.forceY(height / 2).strength(0.05))
                .force("center", d3.forceCenter(width / 2, height / 2))
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
                .force("x", d3.forceX(width / 2).strength(0.05))
                .force("y", d3.forceY(height / 2).strength(0.05))
                .force("center", d3.forceCenter(width / 2, height / 2))
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

            drawChart(chartData.value, {
                width: chartWidth,
                height: chartWidth,
                decade: chartDecade.value
            })

        }

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
    .active {
        background-color: grey;
        color: white;
    }
</style>