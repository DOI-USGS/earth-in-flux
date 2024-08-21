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
        </template>
        <template #figures>
            <div id="chart-container" class="maxWidth" ref="chart"></div>
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
    const chart = ref(null);
    let chartHeight = window.innerHeight * 0.8;
    let bubbleChartDimensions;
    const bubbleChartTitle = 'Title of chart';
    let bubbleChartBounds;
    let barChartDimensions;
    const barChartTitle = 'Title of chart';
    let barChartBounds;
    let timelineChartDimensions;
    const timelineChartTitle = 'Title of chart';
    let timelineChartBounds;
    let chartDecades = ref(null);
    const chartDecade = ref(null);
    const simulation = ref(null);
    const defaultRectOpacity = 0.5;
    const bodyCSS = window.getComputedStyle(document.body);
    const bkgdColor = bodyCSS.getPropertyValue('--color-background');    
    const defaultGrey = '#CECECE';
    const highlightFillGrey = '#969696';
    const highlightStrokeGrey = '#6c6c6c';

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
                const bubbleChartHeight = chartHeight * 0.7;
                initBubbleChart({
                    width: chart.value.offsetWidth, // outer width, in pixels
                    height: bubbleChartHeight
                });
                const timelineChartHeight = chartHeight * 0.075
                initTimelineChart({
                    width: chart.value.offsetWidth, // outer width, in pixels
                    height: timelineChartHeight,
                    marginLeft: 60, // left margin, in pixels
                    marginTop: 0,
                    marginBottom: 25
                });
                const barChartHeight = chartHeight - bubbleChartHeight - timelineChartHeight
                initBarChart({
                    width: chart.value.offsetWidth, // outer width, in pixels
                    height: barChartHeight,
                    marginLeft: 60, // left margin, in pixels
                    marginTop: 10,
                    marginBottom: 20
                })

                // build chart
                drawBubbleChart(data.value, {
                    decade: chartDecade.value
                })
                drawTimelineChart(data.value, {
                    decade: chartDecade.value
                })
                drawBarChart(data.value, {
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
        bubbleChartDimensions = {
            width,
            height,
            margin: {
                top: marginTop,
                right: marginRight,
                bottom: marginBottom,
                left: marginLeft
            }
        }
        bubbleChartDimensions.boundedWidth = bubbleChartDimensions.width - bubbleChartDimensions.margin.left - bubbleChartDimensions.margin.right
        bubbleChartDimensions.boundedHeight = bubbleChartDimensions.height - bubbleChartDimensions.margin.top - bubbleChartDimensions.margin.bottom

        // draw canvas for chart
        const chartSVG = d3.select("#chart-container")
            .append("svg")
                .attr("viewBox", [0, 0, (bubbleChartDimensions.width), (bubbleChartDimensions.height)].join(' '))
                .attr("width", "100%")
                .attr("height", "100%")
                .attr("id", "bubble-chart-svg")

        // assign role for accessibility
        chartSVG.attr("role", "figure")
            .append("title")
            .text(bubbleChartTitle)

        // Add group for bounds
        bubbleChartBounds = chartSVG.append("g")
            .attr("id", "bubble-chart-bounds")
            .style("transform", `translate(${
                bubbleChartDimensions.margin.left
            }px, ${
                bubbleChartDimensions.margin.top
            }px)`)

        bubbleChartBounds.append("g")
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
            .range([4, chart.value.offsetWidth / 7]);

        // filter data to current decade and filter out decades where species is not present
        data = data.filter(d => d.decade === decade && d.pct_abundance > 0);
        
        const nodes = data.map((d) => ({
            ...d,
            radius: sizeScale(parseFloat(d.pct_abundance)),
            x: bubbleChartDimensions.boundedWidth / 2,
            y: bubbleChartDimensions.boundedHeight / 2
        }));
       
        // set up nodes
        let nodeGroups = bubbleChartBounds.selectAll('.nodes')
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
                // .force("x", d3.forceX(bubbleChartDimensions.boundedWidth / 2).strength(0.05))
                // .force("y", d3.forceY(bubbleChartDimensions.boundedHeight / 2).strength(0.05))
                // .force("center", d3.forceCenter(bubbleChartDimensions.boundedWidth / 2, bubbleChartDimensions.boundedHeight / 2))
                // .force(
                //     "collide",
                //     d3.forceCollide()
                //         .radius((d) => d.radius + 2)
                //         .iterations(1)
                // )
                // .force('charge', d3.forceManyBody().strength(0.01))
                // .alphaMin(0.01)
                // .alpha(0)
                // .alphaDecay(0.005)
                .velocityDecay(0.9)
                .on("tick", ticked);
        } else {
            simulation.value = d3.forceSimulation();
            simulation.value
                .nodes(nodes)
                .force("x", d3.forceX(bubbleChartDimensions.boundedWidth / 2).strength(0.05))
                .force("y", d3.forceY(bubbleChartDimensions.boundedHeight / 2).strength(0.05))
                // .force("center", d3.forceCenter(bubbleChartDimensions.boundedWidth / 2, bubbleChartDimensions.boundedHeight / 2))
                .force(
                    "collide",
                    d3.forceCollide()
                        .radius((d) => d.radius + 2)
                        .iterations(1)
                )
                .force('charge', d3.forceManyBody().strength(0.01))
                // .alphaMin(0.01)
                // .alphaDecay(0.005)
                .velocityDecay(0.9)
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

    function initTimelineChart({
        width = 640, // outer width, in pixels
        height = width, // outer height, in pixels
        margin = 1, // default margins
        marginTop = margin, // top margin, in pixels
        marginBottom = margin, // left margin, in pixels
        marginLeft = margin, // left margin, in pixels
        marginRight = margin // right margin, in pixels
    }) {
        // set up global chart dimensions
        timelineChartDimensions = {
            width,
            height,
            margin: {
                top: marginTop,
                right: marginRight,
                bottom: marginBottom,
                left: marginLeft
            }
        }
        timelineChartDimensions.boundedWidth = timelineChartDimensions.width - timelineChartDimensions.margin.left - timelineChartDimensions.margin.right
        timelineChartDimensions.boundedHeight = timelineChartDimensions.height - timelineChartDimensions.margin.top - timelineChartDimensions.margin.bottom

        // draw canvas for chart
        const chartSVG = d3.select("#chart-container")
            .append("svg")
                .attr("viewBox", [0, 0, (timelineChartDimensions.width), (timelineChartDimensions.height)].join(' '))
                .attr("width", "100%")
                .attr("height", "100%")
                .attr("id", "timeline-chart-svg")

        // assign role for accessibility
        chartSVG.attr("role", "figure")
            .append("title")
            .text(timelineChartTitle)

        // Add group for bounds
        timelineChartBounds = chartSVG.append("g")
            .attr("id", "timeline-chart-bounds")
            .style("transform", `translate(${
                timelineChartDimensions.margin.left
            }px, ${
                timelineChartDimensions.margin.top
            }px)`)
    }

    function drawTimelineChart(data, {
        decade = 200
    }) {

        // Set up x scale
        const xScale = d3.scaleBand()
            .domain(chartDecades.value)
            .range([0, timelineChartDimensions.boundedWidth])
            .padding(0);

        // add x axis
        const xAxis = timelineChartBounds.append('g')
            .attr("transform", `translate(0,${timelineChartDimensions.boundedHeight})`)
            .call(d3.axisBottom(xScale).tickSize(0))
            // .select(".domain").remove() // remove axis line;

        xAxis.selectAll('path')
            .attr("stroke", bkgdColor);

        xAxis.selectAll('text')
            .attr("class", d => 'axis-text x-axis label' + d)
            .style('font-weight', d => d === decade ? '900' : '200');

        // Add x axis title
        const xAxisLabelYPosition = xAxis.select("text").attr('y')
        const xAxisLabelDy = xAxis.select("text").attr('dy')
        xAxis.append("text")
            .attr("class", "x-axis axis-title")
            .attr("x", 0)
            .attr("y", xAxisLabelYPosition)
            .attr("dy", xAxisLabelDy)
            .style("text-anchor", "end")
            .text('Year')
       
        // draw timeline
        const lineGroup = timelineChartBounds.append("g")
            .attr("id", "timeline-line-group")

        lineGroup.selectAll('line')
            .data(chartDecades.value[0])
            .enter()
            .append("line")
                .attr("class", "timeline")
                .attr("x1",  d => xScale(d) + xScale.bandwidth() / 2)
                .attr("x2", xScale(chartDecades.value.at(-1)) + xScale.bandwidth() / 2)
                .attr("y1", timelineChartDimensions.boundedHeight / 2)
                .attr("y2", timelineChartDimensions.boundedHeight / 2)
                .style("stroke", defaultGrey)
                .style("stroke-width", 1)

        // draw timeline points
        const pointGroup = timelineChartBounds.append("g")
            .attr("id", "timeline-points-group")

        // Add point for each decade
        pointGroup.selectAll('point')
            .data(chartDecades.value)
            .enter()
            .append("circle")
                .attr("class", d => 'point timeline' + d)
                .attr("cx",  d => xScale(d) + xScale.bandwidth() / 2)
                .attr("cy", timelineChartDimensions.boundedHeight / 2)
                .attr("r", 8)
                .style('fill', d => d === decade ? highlightFillGrey : defaultGrey)
                .style("stroke", d => d === decade ? highlightStrokeGrey : defaultGrey)
                .style("stroke-width", 2)

        // add overlays for chart
        const overlayTimelineGroup = timelineChartBounds.append("g")
            .attr("id", "overlay-timeline-group")

        // Add overlay rectangles over chart for interaction
        overlayTimelineGroup.selectAll('overlays')
            .data(chartDecades.value)
            .enter()
            .append('rect')
                .attr('class', d => 'overlay decade' + d)
                .attr('x', d => xScale(d))
                .attr('y', - timelineChartDimensions.margin.top)
                .attr('height', timelineChartDimensions.boundedHeight + timelineChartDimensions.margin.top)
                .attr('width', xScale.bandwidth())
                .style('fill', bkgdColor)
                .style("opacity", d => d === decade ? 0 : defaultRectOpacity)
                .style("stroke", bkgdColor)
                .style("stroke-opacity", d => d === decade ? 0 : defaultRectOpacity)
                .style("stroke-width", 0.5)
                .on('mouseover', (event, d) => {
                    mouseoverTimelineBar(d)
                })

        // add overlays for axis
        const overlayAxisGroup = timelineChartBounds.append("g")
            .attr("id", "overlay-axis-group")

        // Add overlay rectangles over axis for interaction
        overlayAxisGroup.selectAll('overlays')
            .data(chartDecades.value)
            .enter()
            .append('rect')
                .attr('class', d => 'overlay decade' + d)
                .attr('x', d => xScale(d))
                .attr('y', timelineChartDimensions.boundedHeight)
                .attr('height', timelineChartDimensions.margin.bottom)
                .attr('width', xScale.bandwidth())
                .style('fill', 'transparent')
                .style("stroke-width", 0.5)
                .on('mouseover', (event, d) => {
                    mouseoverTimelineBar(d)
                })
    }

    function initBarChart({
        width = 640, // outer width, in pixels
        height = width, // outer height, in pixels
        margin = 1, // default margins
        marginTop = margin, // top margin, in pixels
        marginBottom = margin, // left margin, in pixels
        marginLeft = margin, // left margin, in pixels
        marginRight = margin // right margin, in pixels
    }) {
        // set up global chart dimensions
        barChartDimensions = {
            width,
            height,
            margin: {
                top: marginTop,
                right: marginRight,
                bottom: marginBottom,
                left: marginLeft
            }
        }
        barChartDimensions.boundedWidth = barChartDimensions.width - barChartDimensions.margin.left - barChartDimensions.margin.right
        barChartDimensions.boundedHeight = barChartDimensions.height - barChartDimensions.margin.top - barChartDimensions.margin.bottom

        // draw canvas for chart
        const chartSVG = d3.select("#chart-container")
            .append("svg")
                .attr("viewBox", [0, 0, (barChartDimensions.width), (barChartDimensions.height)].join(' '))
                .attr("width", "100%")
                .attr("height", "100%")
                .attr("id", "bar-chart-svg")

        // assign role for accessibility
        chartSVG.attr("role", "figure")
            .append("title")
            .text(barChartTitle)

        // Add group for bounds
        barChartBounds = chartSVG.append("g")
            .attr("id", "bar-chart-bounds")
            .style("transform", `translate(${
                barChartDimensions.margin.left
            }px, ${
                barChartDimensions.margin.top
            }px)`)
    }

    function drawBarChart(data, {
        decade = 200
    }) {
        // sort data by bar order, so species with shared color plot together
        data.sort((a,b) => (a.bar_order > b.bar_order) ? 1 : ((b.bar_order > a.bar_order) ? -1 : 0))
        

        // get unique species
        const chartSpecies = d3.union(d3.map(data, d => d.species_id));

        // stack data for rectangles
        const stackedData = d3.stack()
            .keys(chartSpecies)
            .value(([, D], key) => D.get(key)['pct_abundance']) // get value for each series key and stack
            .order(d3.stackOrderNone)
            (d3.index(data, d => d.decade, d => d.species_id));

        // Set up x scale
        const xScale = d3.scaleBand()
            .domain(chartDecades.value)
            .range([0, barChartDimensions.boundedWidth])
            .padding(0);

        // // add x axis
        // const xAxis = barChartBounds.append('g')
        //     .attr("transform", `translate(0,${barChartDimensions.boundedHeight})`)
        //     .call(d3.axisBottom(xScale).tickSize(0));

        // xAxis.selectAll('text')
        //     .attr("class", d => 'axis-text x-axis label' + d)
        //     .style('font-weight', d => d === decade ? '900' : '200');

        // Set up y scale
        const yScale = d3.scaleLinear()
            .domain([0, d3.max(stackedData, d => d3.max(d, d => d[1]))])
            .range([barChartDimensions.boundedHeight, 0]);

        // add y axis
        const yAxis = barChartBounds.append('g')
            .call(d3.axisLeft(yScale)
                .ticks(5)
                .tickFormat(d => d + '%'));

        yAxis.selectAll('text')
            .attr("class", 'axis-text y-axis');

        // set up color scale
        let speciesColors = []
        chartSpecies.forEach(species => {
            const filteredData = data.filter(d => d.species_id === species)
            speciesColors.push(filteredData[0].hexcode)
        })
        const colorScale = d3.scaleOrdinal()
            .domain(chartSpecies)
            .range(speciesColors);
        
        // draw chart
        const rectGroup = barChartBounds.append("g")
            .attr("id", "bar-group")
        // Add subgroup for each category of data
        const speciesRectGroups = rectGroup.selectAll('rects')
            .data(stackedData, d => d.key)
            .enter()
            .append('g')
            .attr("id", d => d.key)

        // Add rectangles for each decade to each species group
        speciesRectGroups.selectAll('rect')
            .data(D => D.map(d => (d.key = D.key, d)))
            .enter().append('rect')
                .attr("class", d => d.key + ' ' + d.data[0].replace(" ", "_"))
                .attr('x', d => xScale(d.data[0]))
                .attr('y', d => yScale(d[1]))
                .attr('height', d => yScale(d[0]) - yScale(d[1]))
                .attr('width', xScale.bandwidth())
                .style("fill", d => colorScale(d.key))
                .style("stroke", d => colorScale(d.key))
                .style("stroke-width", 0.5);

        const overlayGroup = barChartBounds.append("g")
            .attr("id", "overlay-group")

        overlayGroup.selectAll('overlays')
            .data(chartDecades.value)
            .enter()
            .append('rect')
                .attr('class', d => 'overlay decade' + d)
                .attr('x', d => xScale(d))
                .attr('y', - barChartDimensions.margin.top)
                .attr('height', barChartDimensions.boundedHeight + barChartDimensions.margin.top)
                .attr('width', xScale.bandwidth())
                .style('fill', bkgdColor)
                .style("opacity", d => d === decade ? 0 : defaultRectOpacity)
                .style("stroke", bkgdColor)
                .style("stroke-opacity", d => d === decade ? 0 : defaultRectOpacity)
                .style("stroke-width", 0.5)
                .on('mouseover', (event, d) => {
                    mouseoverTimelineBar(d)
                })
    }

    function mouseoverTimelineBar(decade) {
        // If mousedover decade does not match current chartDecade
        if (decade !== chartDecade.value) {
            // Update chart decade
            chartDecade.value = decade;
            // set default styling for points and overlay rectangles
            d3.selectAll('.point')
                .style('fill', defaultGrey)
                .style('stroke', defaultGrey)
            d3.selectAll('.overlay')
                .style('opacity', defaultRectOpacity)
                .style('stroke-opacity', defaultRectOpacity)
            d3.selectAll('.axis-text.x-axis')
                .style('font-weight', '200')
            // set styling for specific decade
            d3.selectAll('.timeline' + decade)
                .style('fill', highlightFillGrey)
                .style('stroke', highlightStrokeGrey)
            d3.selectAll('.decade' + decade)
                .style('opacity', 0)
                .style('stroke-opacity', 0)
            d3.selectAll('.label' + decade)
                .style('font-weight', '900')
            // update bubble chart
            drawBubbleChart(data.value, {decade: chartDecade.value})
        }
    }
</script>

<style lang="scss">
    .axis-text {
        font-size: 1.8rem;
        font-family: var(--default-font);
        user-select: none;
    }
    .axis-title {
        font-size: 1.8rem;
        font-family: var(--default-font);
        font-weight: 900;
        fill: var(--color-text);
        user-select: none;
    }
</style>