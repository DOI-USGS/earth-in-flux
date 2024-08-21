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
    import { onMounted, ref, reactive } from "vue";
    import * as d3 from 'd3';
    import VizSection from '@/components/VizSection.vue';

    // define props
    defineProps({
        text: { type: Object }
    })

    // global variables
    const publicPath = import.meta.env.BASE_URL;
    const dataFile = 'fish_as_food_climate.csv' //'fish_as_food_climate_test.csv'
    const data = ref();
    const chart = ref(null);
    // const families = ref();
    let chartDimensions;
    const chartTitle = 'Title of chart';
    let chartSVG;
    let chartBounds;
    let xScale;
    let xAxis;
    let yScale;
    let yAxis;
    let widthScale;
    let colorScale;
    // Create a reactive object to track expanded families
    // const expandedFamilies = reactive({});


    // Behavior on mounted (functions called here)
    // Load data and then make chart
    onMounted(async () => {
        try {
            await loadDatasets();
            if (data.value.length > 0) {
                // families.value = Array.from(new Set(data.value.map(d => d.family)));
                // // Initialize the expandedFamilies object
                // families.value.forEach(family => {
                //     expandedFamilies[family] = false;
                // });
                // expandedFamilies["A"] = !expandedFamilies["A"]
                // console.log(expandedFamilies)
                initChart({
                    width: chart.value.offsetWidth,
                    height: 1600,
                    margin: 20,
                    marginBottom: 30,
                    marginLeft: 200});
                drawChart(data.value);
            } else {
                console.error('Error loading data');
            }
        } catch (error) {
            console.error('Error during component mounting', error);
        }
    });

    async function loadDatasets() {
        try {
            data.value = await loadData(dataFile);
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

    function initChart({
        width = 500, // outer width, in pixels
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
        chartSVG = d3.select("#chart-container")
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

        // Initialize scales
        initXScale()
        initYScale()

        // Inititalize gradients
        initGradients()
    }

    function initXScale() {
        // scale for the x axis (domain set in `drawChart()`)
        xScale = d3.scaleLinear()
            .range([0, chartDimensions.boundedWidth]);

        // add group for x axis
        xAxis = chartBounds.append("g")
            .attr("id", "x-axis")
            .attr("class", "axis")
            .attr("transform", `translate(0,${chartDimensions.boundedHeight})`)
            .attr("aria-hidden", true) // hide from screen reader

        // generate x axis
        // xAxis
        //     .call(d3.axisBottom(xScale).tickSize(0).tickPadding(10))
            // .select(".domain").remove() // remove axis line
        
        // add placeholder for x axis title (title text set in drawChart())
        xAxis
            .append("text")
            .attr("class", "x-axis axis-title")
            .attr("x", -chartDimensions.boundedWidth / 2)
            .style("text-anchor", "middle")
            .attr("role", "presentation")
            .attr("aria-hidden", true)
    }

    function initYScale() {
        // scale for y axis (domain set in `drawChart()`)
        yScale = d3.scaleBand()
            .range([chartDimensions.boundedHeight, 0])
            .padding(0.1)

        // add group for y axis
        yAxis = chartBounds.append("g")
            .attr("id", "y-axis")
            .attr("class", "axis")
            .attr("aria-hidden", true)

        // generate y axis
        // yAxis
        //     .call(d3.axisLeft(yScale).tickSize(2))
            // .select(".domain").remove() // remove axis line

        // add placeholder for y axis title (title text set in drawChart())
        yAxis
            .append("text")
            .attr("class", "y-axis axis-title")
            .attr("x", -chartDimensions.boundedHeight / 2)
            .attr("transform", "rotate(-90)")
            .style("text-anchor", "middle")
            .attr("role", "presentation")
            .attr("aria-hidden", true)
    }

    function initGradients() {
        let colors = {warm: '#FF7256', cool: '#5CACEE', cold: '#36648B'}
        Object.keys(colors).forEach(color => {
            let lg_decreasing = chartSVG.append("defs").append("linearGradient")
                .attr("id", color + "_gradient_decreasing")
                .attr("x1", "0%")
                .attr("x2", "100%")
                .attr("y1", "0%")
                .attr("y2", "0%");
            lg_decreasing.append("stop")
                .attr("offset", "40%")
                .style("stop-color", colors[color])
                .style("stop-opacity", 1)

            lg_decreasing.append("stop")
                .attr("offset", "100%")
                .style("stop-color", colors[color])
                .style("stop-opacity", 0)

            let lg_increasing = chartSVG.append("defs").append("linearGradient")
                .attr("id", color + "_gradient_increasing")
                .attr("x1", "100%")
                .attr("x2", "0%")
                .attr("y1", "0%")
                .attr("y2", "0%");
            lg_increasing.append("stop")
                .attr("offset", "40%")
                .style("stop-color", colors[color])
                .style("stop-opacity", 1)

            lg_increasing.append("stop")
                .attr("offset", "100%")
                .style("stop-color", colors[color])
                .style("stop-opacity", 0)
        })
    }

    function initWidthScale(low, high) {
        widthScale = d3.scaleLinear()
            .domain([0, 1])
            .range([low, high]);
    }

    function initColorScale(data) {
        colorScale = d3.scaleOrdinal()
            .domain(data)
            .range(data.map(item => getColor(item)));
    }

    function getColor(item) {
        let itemColor;
        switch(item) {
            case 'warm':
                itemColor = '#FF7256';
                break;
            case 'cool':
                itemColor = '#5CACEE';
                break;
            case 'cold':
                itemColor = '#36648B';
                break;
            default:
                itemColor =  '#B3B3B3';
        }
        return itemColor;
    }

    function drawChart(data) {
        // accessor functions
        const yAccessor = d => d.species
        // const yAccessor_family = d => d.family
        const xAccessor = d => d.cvi
        // const xAccessor_family = d => d.cvi_family
        const x0Accessor = d => d.cvi_2030
        const x1Accessor = d => d.cvi_2075
        const x0Accessor_family = d => d.cvi_2030_family
        const x1Accessor_family = d => d.cvi_2075_family
        const widthAccessor = d => d.position
        const colorAccessor = d => d.thermal_guild
        const identifierAccessor = d => d.family + '_' + d.species

        // to get dynamic
        // need key for data
        // need enter update exit pattern with transitions
        // need to transition chart height, yscale, yaxis


        // set domain for xScale
        xScale
            .domain([0, 0.3]) // d3.max([d3.max(data, x0Accessor), d3.max(data, x1Accessor)])])
        xAxis
            .call(d3.axisBottom(xScale).tickSize(0).tickPadding(10))

        xAxis
            .selectAll("text")
            .attr("class", "axis-text")
        // set domain for yScale
        // let yDomain = []
        // families.value.map(family => expandedFamilies[family] ? yDomain.push(data.filter(d => d.family === family).map(d => d.species)) : yDomain.push(family))
        // yDomain = yDomain.flat()
        // yScale
        //     .domain(yDomain)
        yScale
            .domain(d3.union(data.map(yAccessor))) //.sort(d3.ascending)))
        yAxis
            .call(d3.axisLeft(yScale).tickSize(2))
        
        yAxis
            .selectAll("text")
            .attr("class", "axis-text")

        // set up width scale
        const radiusPosition0 = 1
        const strokeRatio = 0.3
        const strokeWidth1 = strokeRatio * yScale.bandwidth() / 2
        const radiusPosition1 = (1 - strokeRatio) * yScale.bandwidth() / 2
        initWidthScale(radiusPosition0, radiusPosition1)

        // set up area function
        const area = d3.area()
            .x(d => xScale(xAccessor(d)))
            .y0(d => yScale(yAccessor(d)) - widthScale(widthAccessor(d)))
            .y1(d => yScale(yAccessor(d)) + widthScale(widthAccessor(d)));
            // .x(d => expandedFamilies[d.family] ? xScale(xAccessor(d)) : xScale(xAccessor_family(d)))
            // .y0(d => expandedFamilies[d.family] ? yScale(yAccessor(d)) - widthScale(widthAccessor(d)) : yScale(yAccessor_family(d)) - widthScale(widthAccessor(d)))
            // .y1(d => expandedFamilies[d.family] ? yScale(yAccessor(d)) + widthScale(widthAccessor(d)) : yScale(yAccessor_family(d)) + widthScale(widthAccessor(d)));

        // set up color scale
        const colorCategories = Array.from(new Set(data.map(colorAccessor)))
        initColorScale(colorCategories)

        // set up area data
        const areaCategories = Array.from(new Set(data.map(yAccessor)))
        const areaData = areaCategories.map(areaCategory => {
            return [
                {
                    species: areaCategory,
                    family: data.filter(d => d.species === areaCategory)[0].family,
                    cvi: x0Accessor(data.filter(d => d.species === areaCategory)[0]),
                    cvi_decreasing: x0Accessor(data.filter(d => d.species === areaCategory)[0]) > x1Accessor(data.filter(d => d.species === areaCategory)[0]),
                    cvi_family: x0Accessor_family(data.filter(d => d.species === areaCategory)[0]),
                    position: 0,
                    thermal_guild: colorAccessor(data.filter(d => d.species === areaCategory)[0])
                },
                {
                    species: areaCategory,
                    family: data.filter(d => d.species === areaCategory)[0].family,
                    cvi: x1Accessor(data.filter(d => d.species === areaCategory)[0]),
                    cvi_decreasing: x0Accessor(data.filter(d => d.species === areaCategory)[0]) > x1Accessor(data.filter(d => d.species === areaCategory)[0]),
                    cvi_family: x1Accessor_family(data.filter(d => d.species === areaCategory)[0]),
                    position: 1,
                    thermal_guild: colorAccessor(data.filter(d => d.species === areaCategory)[0])
                }
            ]
        })
            
        // draw chart
        chartBounds.append("g")
            .attr("id", "areas")
            .selectAll('.area')
                .data(areaData)
                .enter()
                .append('path')
                    .attr("id", d => 'area-2030-' + identifierAccessor(d))
                    .attr('class', "area")
                    .attr('d', d => area(d))
                    .attr('fill', d => d[0].cvi_decreasing ? `url(#${d[0].thermal_guild}_gradient_decreasing)`: `url(#${d[0].thermal_guild}_gradient_increasing)`)//d => colorScale(colorAccessor(d[0])))
                    .style("opacity", 1)
                    // .on("click", function(event, d) {
                    //     const clickedFamily = d[0].family
                    //     expandedFamilies[clickedFamily] = !expandedFamilies[clickedFamily]
                    //     let yDomain = []
                    //     families.value.map(family => expandedFamilies[family] ? yDomain.push(data.filter(d => d.family === family).map(d => d.species)) : yDomain.push(family))
                    //     yDomain = yDomain.flat()
                    //     const yBandwidth = yScale.bandwidth()
                    //     chartDimensions.boundedHeight = yBandwidth * yDomain.length;
                    //     chartDimensions.height = chartDimensions.boundedHeight + chartDimensions.margin.top + chartDimensions.margin.bottom;
                    //     chartSVG.attr("viewBox", [0, 0, (chartDimensions.width), (chartDimensions.height)].join(' '))
                    //     yScale.range([chartDimensions.boundedHeight * 2, 0])
                    //     yAxis
                    //         .call(d3.axisLeft(yScale).tickSize(2))
                    //     drawChart(data)
                    // });

        chartBounds.append("g")
            .attr("id", "points-2030")
            .attr("class", "points points_2030")    
            .selectAll('points')
                .data(data)
                .enter()
                .append("circle")
                    .attr("id", d => 'point-2030-' + identifierAccessor(d))
                    .attr("class", "point")
                    .attr("cx", d => xScale(x0Accessor(d)))
                    .attr("cy", d => yScale(yAccessor(d)))
                    .attr("r", radiusPosition0)
                    .style("stroke", d => colorScale(colorAccessor(d)))
                    .style("fill", "white")

        chartBounds.append("g")
            .attr("id", "points-2075")
            .attr("class", "points points_2075")    
            .selectAll('points')
                .data(data)
                .enter()
                .append("circle")
                    .attr("id", d => 'point-2075-' + identifierAccessor(d))
                    .attr("class", "point")
                    .attr("cx", d => xScale(x1Accessor(d)))
                    .attr("cy", d => yScale(yAccessor(d)))
                    .attr("r", radiusPosition1 - strokeWidth1 / 2)
                    .style("stroke", d => colorScale(colorAccessor(d)))
                    .style("stroke-width", strokeWidth1)
                    .style("fill", "white")

        // chartBounds.append("g")
        //     .attr("id", "eyes-2075")
        //     .attr("class", "eyes eyes_2075")    
        //     .selectAll('eyes')
        //         .data(data)
        //         .enter()
        //         .append("circle")
        //             .attr("id", d => 'point-2075-' + identifierAccessor(d))
        //             .attr("class", "point")
        //             .attr("cx", d => {
        //                 return x1Accessor(d) > x0Accessor(d) ? (xScale(x1Accessor(d)) + radiusPosition1 / 4) : (xScale(x1Accessor(d)) - radiusPosition1 / 4)
        //             })
        //             .attr("cy", d => yScale(yAccessor(d)) + radiusPosition1 / 4)
        //             .attr("r", radiusPosition1 / 2)
        //             .style("stroke", 'none')
        //             .style("fill", "black")

    }
</script>

<style lang="scss">
    .axis-text {
        font-size: 1.8rem;
        font-family: var(--default-font);
        user-select: none;
    }
</style>