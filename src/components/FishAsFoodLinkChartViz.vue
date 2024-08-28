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
            <button
                  aria-pressed="!scalePercent"
                  class="button"
                  :text="scaleType"
                  @click="toggleScale"
                >
                  {{ scaleType }}
            </button>
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
    import { computed, onMounted, reactive, ref } from "vue"; //, reactive
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
    const families = ref();
    const chart = ref(null);
    const scalePercent = ref(false);
    let chartDimensions;
    const chartTitle = 'Title of chart';
    let chartSVG;
    let chartBounds;
    let xScale;
    let xAxisBottom;
    let yScale;
    let yAxis;
    let widthScale;
    let colorScale;
    // Create a reactive object to track expanded families
    const expandedFamilies = reactive({});

    // set up filtered chart data as computed property
    const scaleType = computed(() => {
        return scalePercent.value ? 'Percent change' : 'Change'
    });


    // Behavior on mounted (functions called here)
    // Load data and then make chart
    onMounted(async () => {
        try {
            await loadDatasets();
            if (data.value.length > 0) {
                families.value = Array.from(new Set(data.value.map(d => d.family)));
                // Initialize the expandedFamilies object
                families.value.forEach(family => {
                    expandedFamilies[family] = false;
                });
                // expandedFamilies["Cyprinidae"] = !expandedFamilies["Cyprinidae"]
                expandedFamilies["Salmonidae"] = !expandedFamilies["Salmonidae"]

                initChart({
                    width: chart.value.offsetWidth,
                    height: 1800,
                    margin: 20,
                    marginBottom: 80,
                    marginLeft: 200});

                drawChart(data.value, scalePercent.value);
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
            const data = await d3.csv(publicPath + fileName, d => {
                d.cvi_2030 = +d.cvi_2030;
                d.cvi_2075 = +d.cvi_2075;
                d.cvi_2030_family = +d.cvi_2030_family;
                d.cvi_2075_family = +d.cvi_2075_family;
                return d;
            });
            return data;
        } catch (error) {
            console.error(`Error loading data from ${fileName}`, error);
            return [];
        }
    }

    function toggleScale() {
        scalePercent.value = !scalePercent.value
        drawChart(data.value, scalePercent.value)
    }

    function initChart({
        width = 500, // outer width, in pixels
        // height = width, // outer height, in pixels
        margin = 1, // default margins
        marginTop = margin, // top margin, in pixels
        marginBottom = margin, // left margin, in pixels
        marginLeft = margin, // left margin, in pixels
        marginRight = margin // right margin, in pixels
    }) {
        // set up global chart dimensions
        chartDimensions = {
            width,
            //height,
            margin: {
                top: marginTop,
                right: marginRight,
                bottom: marginBottom,
                left: marginLeft
            }
        }
        chartDimensions.boundedWidth = chartDimensions.width - chartDimensions.margin.left - chartDimensions.margin.right
        //chartDimensions.boundedHeight = chartDimensions.height - chartDimensions.margin.top - chartDimensions.margin.bottom

        // draw canvas for chart
        chartSVG = d3.select("#chart-container")
            .append("svg")
                //.attr("viewBox", [0, 0, (chartDimensions.width), (chartDimensions.height)].join(' '))
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

        // Initialize axes
        initXAxis()
        initYAxis()

        // Initialize scales
        initXScale()
        initYScale()

        // Inititalize gradients
        initGradients()

        // Add groups for visual elements
        chartBounds.append("g")
            .attr("class", "areas")
        chartBounds.append("g")
            .attr("class", "points_2030")
        chartBounds.append("g")
            .attr("class", "points_2075")
        chartBounds.append("g")
            .attr("class", "rects")
    }

    function initXScale() {
        // scale for the x axis (domain set in `drawChart()`)
        xScale = d3.scaleLinear()
            .range([0, chartDimensions.boundedWidth]);
    }
    function initXAxis() {
        // add group for x axis
        xAxisBottom = chartBounds.append("g")
            .attr("id", "x-axis")
            .attr("class", "axis")
            //.attr("transform", `translate(0,${chartDimensions.boundedHeight})`)
            .attr("aria-hidden", true) // hide from screen reader

        // generate x axis
        // xAxisBottom
        //     .call(d3.axisBottom(xScale).tickSize(0).tickPadding(10))
            // .select(".domain").remove() // remove axis line
        
        // // add placeholder for x axis title (title text set in drawChart())
        // xAxisBottom
        //     .append("text")
        //     .attr("class", "x-axis axis-title")
        //     .attr("x", -chartDimensions.boundedWidth / 2)
        //     .style("text-anchor", "middle")
        //     .attr("role", "presentation")
        //     .attr("aria-hidden", true)
    }

    function initYScale() {
        // scale for y axis (domain set in `drawChart()`)
        yScale = d3.scaleBand()
            // .range([chartDimensions.boundedHeight, 0])
            .padding(0.1)


    }
    function initYAxis() {
        // add group for y axis
        yAxis = chartBounds.append("g")
            .attr("id", "y-axis")
            .attr("class", "axis")
            .attr("aria-hidden", true)

        // generate y axis
        // yAxis
        //     .call(d3.axisLeft(yScale).tickSize(2))
            // .select(".domain").remove() // remove axis line

        // // add placeholder for y axis title (title text set in drawChart())
        // yAxis
        //     .append("text")
        //     .attr("class", "y-axis axis-title")
        //     .attr("x", -chartDimensions.boundedHeight / 2)
        //     .attr("transform", "rotate(-90)")
        //     .style("text-anchor", "middle")
        //     .attr("role", "presentation")
        //     .attr("aria-hidden", true)
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

    function drawChart(data, scalePercent) {
        // accessor functions
        const yAccessor = d => expandedFamilies[d.family] ? d.species : d.family //d.species
        // const yAccessor_family = d => d.family
        const xAccessor = d => d.cvi
        // const xAccessor_family = d => d.cvi_family
        const x0Accessor = d => scalePercent ? 0 : d.cvi_2030
        const x1Accessor = d => scalePercent ? (d.cvi_2075 - d.cvi_2030)/ d.cvi_2030 : d.cvi_2075
        const x0Accessor_family = d => scalePercent ? 0 : d.cvi_2030_family
        const x1Accessor_family = d => scalePercent ? (d.cvi_2075_family - d.cvi_2030_family)/ d.cvi_2030_family : d.cvi_2075_family
        const widthAccessor = d => d.position
        const colorAccessor = d => d.thermal_guild
        const identifierAccessor = d => expandedFamilies[d.family] ? d.family + '_' + d.species.replace(/ /g,"_") : d.family;

        // to get dynamic
        // need key for data
        // need enter update exit pattern with transitions
        // need to transition chart height, yscale, yaxis


        // set domain for xScale
        if (scalePercent) {
            const maxVal = d3.max([Math.abs(d3.max(data, x1Accessor)), Math.abs(d3.min(data, x1Accessor))])
            xScale
                .domain([-maxVal, maxVal])
                .nice()
        } else {
            const maxVal = d3.max([d3.min(data, x0Accessor), d3.max(data, x0Accessor), d3.min(data, x1Accessor), d3.max(data, x1Accessor)])
            xScale
                .domain([0, Math.round(maxVal * 10) / 10])
                .nice()
        }
        
        xAxisBottom.transition(getUpdateTransition())
            .call(d3.axisBottom(xScale).tickSize(0).tickPadding(10));

        xAxisBottom
            .selectAll("text")
            .attr("class", "axis-text")

        const xAxisBottomLabelYPosition = xAxisBottom.select("text").attr('y')
        const xAxisBottomLabelDy = xAxisBottom.select("text").attr('dy')
        xAxisBottom.append("text")
            .attr("class", "x-axis axis-title")
            .attr("x", chartDimensions.boundedWidth / 2)
            .attr("y", xAxisBottomLabelYPosition * 4)
            .attr("dy", xAxisBottomLabelDy)
            .style("text-anchor", "middle")
            .text(scalePercent ? 'Percent change in harvest-weighted climate vulnerability, 2030-2075' : 'Change in harvest-weighted climate vulnerability, 2030-2075')
        xAxisBottom.append("text")
            .attr("class", "x-axis axis-subtitle")
            .attr("x", 0)
            .attr("y", xAxisBottomLabelYPosition * 5.5)
            .attr("dy", xAxisBottomLabelDy)
            .style("text-anchor", "start")
            .text('Less vulnerable')
        xAxisBottom.append("text")
            .attr("class", "x-axis axis-subtitle")
            .attr("x", chartDimensions.boundedWidth)
            .attr("y", xAxisBottomLabelYPosition * 5.5)
            .attr("dy", xAxisBottomLabelDy)
            .style("text-anchor", "end")
            .text('More vulnerable')

        // Remove axix line and labels
        // xAxisBottom.select(".domain").remove()
        // xAxisBottom.call(d3.axisBottom(xScale).tickValues([]))

        // set domain for yScale
        let yDomain = []
        families.value.map(family => expandedFamilies[family] ? yDomain.push(data.filter(d => d.family === family).map(d => d.species)) : yDomain.push(family))
        yDomain = yDomain.flat()
        
        const totalHeight = yDomain.length * 25;
        chartDimensions.boundedHeight = totalHeight
        chartDimensions.height = chartDimensions.boundedHeight + chartDimensions.margin.top + chartDimensions.margin.bottom
        chartSVG
            .transition(getUpdateTransition())
            .attr("viewBox", [0, 0, (chartDimensions.width), (chartDimensions.height)].join(' '))
        
        // Set range and domain for y scale
        yScale
            .range([chartDimensions.boundedHeight, 0])
            .domain(yDomain)
        //yScale
        //    .domain(d3.union(data.map(yAccessor))) //.sort(d3.ascending)))
        
        // set y position for xAxisBottom
        xAxisBottom
            .transition(getUpdateTransition())
            .attr("transform", `translate(0,${chartDimensions.boundedHeight})`)

        yAxis
            .transition(getUpdateTransition())
            .call(d3.axisLeft(yScale).tickSize(2))
        
        yAxis
            .selectAll("text")
            .attr("class", d => expandedFamilies[d] === false ? "axis-text family" : "axis-text species")

        const yAxisLabelXPosition = yAxis.select("text").attr('x')
        const yAxisLabelDx = yAxis.select("text").attr('dx')
        // yAxis.append("text")
        //     .attr("class", "y-axis axis-title")
        //     .attr("y", 0)
        //     .attr("x", yAxisLabelXPosition)
        //     .attr("dx", yAxisLabelDx)
        //     .style("text-anchor", "end")
        //     .text('Species')

        // set up width scale
        const radiusPosition0 = 1
        const strokeRatio = 0.3
        const strokeWidth1 = strokeRatio * yScale.bandwidth() / 2
        const radiusPosition1 = (1 - strokeRatio) * yScale.bandwidth() / 2
        initWidthScale(radiusPosition0, radiusPosition1)

        // set up area function
        const area = d3.area()
            .x(d => xScale(xAccessor(d)))
            .y0(d => yScale(yAccessor(d)) + yScale.bandwidth() / 2 - widthScale(widthAccessor(d)))
            .y1(d => yScale(yAccessor(d)) + yScale.bandwidth() / 2+ widthScale(widthAccessor(d)));
            // .x(d => expandedFamilies[d.family] ? xScale(xAccessor(d)) : xScale(xAccessor_family(d)))
            // .y0(d => expandedFamilies[d.family] ? yScale(yAccessor(d)) - widthScale(widthAccessor(d)) : yScale(yAccessor_family(d)) - widthScale(widthAccessor(d)))
            // .y1(d => expandedFamilies[d.family] ? yScale(yAccessor(d)) + widthScale(widthAccessor(d)) : yScale(yAccessor_family(d)) + widthScale(widthAccessor(d)));

        // set up color scale
        const colorCategories = Array.from(new Set(data.map(colorAccessor)))
        initColorScale(colorCategories)

        // set up area data
        const areaCategories = Array.from(new Set(data.map(d => d.species))) //Array.from(new Set(data.map(yAccessor)))
        const areaData = areaCategories.map(areaCategory => {
            const species = areaCategory //expandedFamilies[areaCategory] === false ? null : areaCategory;
            const family = data.filter(d => d.species === species)[0].family //expandedFamilies[areaCategory] === false ? areaCategory : data.filter(d => d.species === species)[0].family;
            const cvi_2030 = expandedFamilies[family] ? x0Accessor(data.filter(d => d.species === species)[0]) : x0Accessor_family(data.filter(d => d.family === family)[0]) //species ? x0Accessor(data.filter(d => d.species === species)[0]) : x0Accessor_family(data.filter(d => d.family === family)[0])
            const cvi_2075 = expandedFamilies[family] ? x1Accessor(data.filter(d => d.species === species)[0]) : x1Accessor_family(data.filter(d => d.family === family)[0])
            const cvi_decreasing = cvi_2030 > cvi_2075
            const thermal_guild = colorAccessor(data.filter(d => d.species === species)[0]) //species ? colorAccessor(data.filter(d => d.species === species)[0]) : colorAccessor(data.filter(d => d.family === family)[0])

            return [
                {
                    species: species,
                    family: family,
                    cvi: cvi_2030,
                    cvi_decreasing: cvi_decreasing,
                    // cvi_family: x0Accessor_family(data.filter(d => d.species === areaCategory)[0]),
                    position: 0,
                    thermal_guild: thermal_guild
                },
                {
                    species: species,
                    family: family,
                    cvi: cvi_2075,
                    cvi_decreasing: cvi_decreasing,
                    // cvi_family: x1Accessor_family(data.filter(d => d.species === areaCategory)[0]),
                    position: 1,
                    thermal_guild: thermal_guild
                }
            ]
        })
            
        // draw chart
        // Enter-update-exit pattern for areas
        let areaGroups = chartBounds.selectAll(".areas")
            .selectAll(".area")
            .data(areaData, d => d[0].species + d[0].family)
        
        const oldAreaGroups = areaGroups.exit()

        oldAreaGroups.selectAll('path')
            .transition(getExitTransition())
            .style("opacity", 0)

        oldAreaGroups.transition(getExitTransition()).remove()
        
        const newAreaGroups = areaGroups.enter().append("g")
            .attr("class", d => "area " + d[0].species)
            .attr("id", d => 'area-group-' + identifierAccessor(d[0]))

        // append paths
        newAreaGroups.append("path")
            .append('path')
            .attr("id", d => 'area-' + identifierAccessor(d[0]))
            .attr('d', null)
            .attr('fill', d => d[0].cvi_decreasing ? `url(#${d[0].thermal_guild}_gradient_decreasing)`: `url(#${d[0].thermal_guild}_gradient_increasing)`)//d => colorScale(colorAccessor(d[0])))
            .style("opacity", 1)

        // update areaGroups to include new paths
        areaGroups = newAreaGroups.merge(areaGroups)

        const areaPaths = areaGroups.select("path")

        // Update paths based on data values
        areaPaths.transition(getUpdateTransition())
            .attr("id", d => 'area-' + identifierAccessor(d[0]))
            .attr('d', d => area(d))
            .attr('fill', d => d[0].cvi_decreasing ? `url(#${d[0].thermal_guild}_gradient_decreasing)`: `url(#${d[0].thermal_guild}_gradient_increasing)`)//d => colorScale(colorAccessor(d[0])))
            .style("opacity", 1)
        
        // chartBounds.append("g")
        //     .attr("id", "areas")
        //     .selectAll('.area')
        //         .data(areaData, d => d[0].species)
        //         .enter()
        //         .append('path')
        //             .attr("id", d => 'area-' + identifierAccessor(d[0]))
        //             .attr('class', "area")
        //             .attr('d', d => area(d))
        //             .attr('fill', d => d[0].cvi_decreasing ? `url(#${d[0].thermal_guild}_gradient_decreasing)`: `url(#${d[0].thermal_guild}_gradient_increasing)`)//d => colorScale(colorAccessor(d[0])))
        //             .style("opacity", 1)
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
                    //     drawChart(data, scalePercent.value)
                    // });

        // Enter-Update-Exit pattern for 2030 points
        let pointGroups2030 = chartBounds.selectAll('.points_2030')
            .selectAll(".point_2030")
            .data(data, d => d.species)
        
        const oldPointGroups2030 = pointGroups2030.exit()

        oldPointGroups2030.selectAll('circle')
            .transition(getExitTransition())
            .style("opacity", 0)

        oldPointGroups2030.transition(getExitTransition()).remove()
        
        const newPointGroups2030 = pointGroups2030.enter().append("g")
            .attr("class", d => "point_2030 " + d.species)
            .attr("id", d => 'point-2030-group-' + identifierAccessor(d))

        // append points
        newPointGroups2030
            .append("circle")
                .attr("id", d => 'point-2030-' + identifierAccessor(d))
                .attr("class", "point_2030")
                .attr("cx", d => expandedFamilies[d.family] ? xScale(x0Accessor(d)) : xScale(x0Accessor_family(d)))
                .attr("cy", d => yScale(yAccessor(d)) + yScale.bandwidth() / 2)
                .attr("r", radiusPosition0)
                .style("stroke", d => colorScale(colorAccessor(d)))
                .style("fill", "white")

        // update pointGroups2030 to include new points
        pointGroups2030 = newPointGroups2030.merge(pointGroups2030)

        const allPoints2030 = pointGroups2030.select("circle")

        // Update points based on data values
        allPoints2030.transition(getUpdateTransition())
            .attr("id", d => 'point-2030-' + identifierAccessor(d))
            .attr("class", "point_2030")
            .attr("cx", d => expandedFamilies[d.family] ? xScale(x0Accessor(d)) : xScale(x0Accessor_family(d)))
            .attr("cy", d => yScale(yAccessor(d)) + yScale.bandwidth() / 2)
            .attr("r", radiusPosition0)
            .style("stroke", d => colorScale(colorAccessor(d)))
            .style("fill", "white")

        // chartBounds.append("g")
        //     .attr("id", "points-2030")
        //     .attr("class", "points points_2030")    
        //     .selectAll('points')
        //         .data(data)
        //         .enter()
        //         .append("circle")
        //             .attr("id", d => 'point-2030-' + identifierAccessor(d))
        //             .attr("class", "point")
        //             .attr("cx", d => xScale(x0Accessor(d)))
        //             .attr("cy", d => yScale(yAccessor(d)) + yScale.bandwidth() / 2)
        //             .attr("r", radiusPosition0)
        //             .style("stroke", d => colorScale(colorAccessor(d)))
        //             .style("fill", "white")

        // Enter-Update-Exit pattern for 2075 points
        let pointGroups2075 = chartBounds.selectAll('.points_2075')
            .selectAll(".point_2075")
            .data(data, d => d.species)
        
        const oldPointGroups2075 = pointGroups2075.exit()

        oldPointGroups2075.selectAll('circle')
            .transition(getExitTransition())
            .style("opacity", 0)

        oldPointGroups2075.transition(getExitTransition()).remove()
        
        const newPointGroups2075 = pointGroups2075.enter().append("g")
            .attr("class", d => "point_2075 " + d.species)
            .attr("id", d => 'point-2030-group-' + identifierAccessor(d))

        // append points
        newPointGroups2075
            .append("circle")
                .attr("id", d => 'point-2075-' + identifierAccessor(d))
                .attr("class", "point_2075")
                .attr("cx", d => expandedFamilies[d.family] ? xScale(x1Accessor(d)) : xScale(x1Accessor_family(d)))
                .attr("cy", d => yScale(yAccessor(d)) + yScale.bandwidth() / 2)
                .attr("r", radiusPosition1 - strokeWidth1 / 2)
                .style("stroke", d => colorScale(colorAccessor(d)))
                .style("stroke-width", strokeWidth1)
                .style("fill", "white")

        // update pointGroups2075 to include new points
        pointGroups2075 = newPointGroups2075.merge(pointGroups2075)

        const allPoints2075 = pointGroups2075.select("circle")

        // Update points based on data values
        allPoints2075.transition(getUpdateTransition())
            .attr("id", d => 'point-2075-' + identifierAccessor(d))
            .attr("class", "point_2075")
            .attr("cx", d => expandedFamilies[d.family] ? xScale(x1Accessor(d)) : xScale(x1Accessor_family(d)))
            .attr("cy", d => yScale(yAccessor(d)) + yScale.bandwidth() / 2)
            .attr("r", radiusPosition1 - strokeWidth1 / 2)
            .style("stroke", d => colorScale(colorAccessor(d)))
            .style("stroke-width", strokeWidth1)
            .style("fill", "white")

        // chartBounds.append("g")
        //     .attr("id", "points-2075")
        //     .attr("class", "points points_2075")    
        //     .selectAll('points')
        //         .data(data)
        //         .enter()
        //         .append("circle")
        //             .attr("id", d => 'point-2075-' + identifierAccessor(d))
        //             .attr("class", "point")
        //             .attr("cx", d => xScale(x1Accessor(d)))
        //             .attr("cy", d => yScale(yAccessor(d)) + yScale.bandwidth() / 2)
        //             .attr("r", radiusPosition1 - strokeWidth1 / 2)
        //             .style("stroke", d => colorScale(colorAccessor(d)))
        //             .style("stroke-width", strokeWidth1)
        //             .style("fill", "white")

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

        // Enter-Update-Exit pattern for overlay rectangles
        let rectGroups = chartBounds.selectAll('.rects')
            .selectAll(".rect")
            .data(data, d => d.species)
        
        const oldRectGroups = rectGroups.exit()

        oldRectGroups.selectAll('rect')
            .transition(getExitTransition())
            .style("opacity", 0)

        oldRectGroups.transition(getExitTransition()).remove()
        
        const newRectGroups = rectGroups.enter().append("g")
            .attr("class", d => "rect " + d.species)
            .attr("id", d => 'rect-group-' + identifierAccessor(d))

        // append points
        newRectGroups
            .append("rect")
                .attr("id", d => 'rect-' + identifierAccessor(d))
                .attr("class", "rect")
                .attr("x", -chartDimensions.margin.left)
                .attr("y", d => yScale(yAccessor(d)))
                .attr("height", yScale.bandwidth())
                .attr("width", chartDimensions.width)
                .style("fill", "transparent")

        // update pointGroups2075 to include new points
        rectGroups = newRectGroups.merge(rectGroups)

        const allRects = rectGroups.select("rect")

        // Update points based on data values
        allRects.transition(getUpdateTransition())
            .attr("id", d => 'rect-' + identifierAccessor(d))
            .attr("class", "rect")
            .attr("x", -chartDimensions.margin.left)
            .attr("y", d => yScale(yAccessor(d)))
            .attr("height", yScale.bandwidth())
            .attr("width", chartDimensions.width)
            .style("fill", "transparent")
        
        allRects
            .on("click", function(event, d) {
                const clickedFamily = d.family
                expandedFamilies[clickedFamily] = !expandedFamilies[clickedFamily]
                // let yDomain = []
                // families.value.map(family => expandedFamilies[family] ? yDomain.push(data.filter(d => d.family === family).map(d => d.species)) : yDomain.push(family))
                // yDomain = yDomain.flat()
                // const yBandwidth = yScale.bandwidth()
                // chartDimensions.boundedHeight = yBandwidth * yDomain.length;
                // chartDimensions.height = chartDimensions.boundedHeight + chartDimensions.margin.top + chartDimensions.margin.bottom;
                // chartSVG.attr("viewBox", [0, 0, (chartDimensions.width), (chartDimensions.height)].join(' '))
                // yScale.range([chartDimensions.boundedHeight * 2, 0])
                // yAxis
                //     .call(d3.axisLeft(yScale).tickSize(2))
                drawChart(data, scalePercent)
            });
    }

    function getUpdateTransition () {
      return d3.transition()
        .duration(1500)
        .ease(d3.easeCubicInOut)
    }
    function getExitTransition() {
      return d3.transition()
        .duration(1500)
        .ease(d3.easeCubicInOut)
    }
</script>

<style lang="scss">
    .axis-text {
        font-size: 1.8rem;
        font-family: var(--default-font);
        user-select: none;
        
    }
    .species {
        font-style: italic;
    }
    .family {
        font-weight: 900;
    }
    .axis-title {
        font-size: 1.8rem;
        font-family: var(--default-font);
        font-weight: 900;
        fill: var(--color-text);
        user-select: none;
    }
    .axis-subtitle {
        font-size: 1.8rem;
        font-family: var(--default-font);
        font-weight: 400;
        font-style: italic;
        fill: var(--color-text);
        user-select: none;
    }
</style>