<template>
    <!---VizSection-->
    <VizSection
        id="cross-section"
        :figures="true"
        :fig-caption="false"
    >
        <template #heading>
            <h2>
                {{ text.heading }}
            </h2>
        </template>
        <template #figures>
            <div id="wildfire-aerosols-grid-container">
                <button id="aerosol-prev" class="flip-button" @click="currentIndex--" :disabled="isFirstImage">
                    <font-awesome-icon :icon="{ prefix: 'fas', iconName: 'arrow-left' }"  class="fa fa-arrow-left"/>
                </button>
                <button id="aerosol-next" class="flip-button" @click="currentIndex++" :disabled="isLastImage">
                    <font-awesome-icon :icon="{ prefix: 'fas', iconName: 'arrow-right' }"  class="fa fa-arrow-right"/>
                </button>
                <div id="aerosol-text-container" class="text-container">
                    <p v-html="currentText" />
                </div>
                <div id="chart-container" ref="chart"></div>
            </div>
        </template>
        <template #figureCaption>
        </template>
        <template #belowExplanation>
        </template>
    </VizSection>
</template>

<script setup>
    import { computed, onMounted, ref } from "vue";
    import * as d3 from 'd3';
    import VizSection from '@/components/VizSection.vue';

    // define props
    const props = defineProps({
        text: { type: Object }
    })

    // global variables
    const publicPath = import.meta.env.BASE_URL;
    const tileDataFile = 'fii_core4particulates.csv';
    const barDataFile = 'fii_core4sugars.csv';
    const scatterDataFile = 'fii_core4biomass.csv';
    const tileData = ref();
    const barData = ref();
    const scatterData = ref();
    const currentIndex = ref(1);
    const nIndices = 3;
    const chart = ref(null);
    const chartTitle = 'Title of chart';
    const chartHeight = window.innerHeight * 0.8;
    let chartWidth;
    let chartDimensions;
    let chartBounds;
    let tileChartDimensions;
    let tileChartBounds;
    let tileColorScale;
    // let xScale;
    // let xAxis;
    let yScale;
    let yAxis;
    let barChartDimensions;
    let barChartBounds;
    let barXScale;
    let barYScale;
    let barXAxis;
    const barColors = {Mannosan: '#000000', Galactosan: '#989898', Levoglucosan: '#c8c8c8'};
    let barColorScale;
    let scatterChartDimensions;
    let scatterChartBounds;
    let scatterXScale;
    const scatterColors = {grass: '#c49051', hardwood: '#3c475a', softwood: '#729C9D'};
    let scatterColorScale;

    const isFirstImage = computed(() => {
        return currentIndex.value === 1;
    });
    const isLastImage = computed(() => {
        return currentIndex.value === nIndices;
    });
    const currentText = computed(() => {
        const selectionString = 'paragraph' + currentIndex.value
        return props.text[selectionString];
    });

    // Behavior on mounted (functions called here)
    // Load data and then make chart
    onMounted(async () => {
        try {
            await loadDatasets({
                dataFiles: [tileDataFile, barDataFile, scatterDataFile], 
                dataRefs: [tileData, barData, scatterData],
                dataNumericFields: [['depth_cm', 'total'], ['picogram_per_mL'], null]
            });

            if (tileData.value.length > 0 && barData.value.length > 0) {
                // initialize chart elements
                chartWidth = chart.value.offsetWidth;
                initChart({
                    width: chart.value.offsetWidth,
                    height: chartHeight,
                    margin: 10
                })

                const tileChartWidth = chartWidth / 3
                initTileChart({
                    width: tileChartWidth,
                    height: chartHeight,
                    margin: 10,
                    marginBottom: 50,
                    marginLeft: 100});

                const barChartWidth = chartWidth / 3
                initBarChart({
                    width: barChartWidth,
                    height: chartHeight,
                    margin: 10,
                    marginBottom: 50,
                    marginLeft: 80,
                    translateX: tileChartWidth});

                const scatterChartWidth = chartWidth - tileChartWidth - barChartWidth
                initScatterChart({
                    width: scatterChartWidth,
                    height: chartHeight,
                    margin: 10,
                    marginBottom: 50,
                    marginLeft: 60,
                    translateX: tileChartWidth + barChartWidth});

                // draw charts
                drawTileChart(tileData.value);
                drawBarChart(barData.value);
                drawScatterChart(scatterData.value)
            } else {
                console.error('Error loading data');
            }
        } catch (error) {
            console.error('Error during component mounting', error);
        }
    });

    async function loadDatasets({dataFiles, dataRefs, dataNumericFields}) {
        try {
            for (let i = 0; i < Math.min(dataFiles.length, dataRefs.length, dataNumericFields.length); i++) {
                dataRefs[i].value = await loadData(dataFiles[i], dataNumericFields[i]);
                console.log(`${dataFiles[i]} data in`);
            }
        } catch (error) {
            console.error('Error loading datasets', error);
        }
    }

    async function loadData(dataFile, dataNumericFields) {
        try {
            const data = await d3.csv(publicPath + dataFile, d => {
                if (dataNumericFields) {
                    dataNumericFields.forEach(numericField => {
                        d[numericField] = +d[numericField]
                    });
                }
                return d;
            });
            return data;
        } catch (error) {
            console.error(`Error loading data from ${dataFile}`, error);
            return [];
        }
    }

    function initChart({
        width = 500, // outer width, in pixels
        height = 500, // outer height, in pixels
        margin = 1, // default margins
        marginTop = margin, // top margin, in pixels
        marginBottom = margin, // left margin, in pixels
        marginLeft = margin, // left margin, in pixels
        marginRight = margin // right margin, in pixels
    }) {
        // set up global chart dimensions, including bounded dimensions
        chartDimensions = {
            width,
            height,
            margin: {
                top: marginTop,
                right: marginRight,
                bottom: marginBottom,
                left: marginLeft
            },
            boundedWidth: width - marginLeft - marginRight,
            boundedHeight: height - marginTop - marginBottom
        }

        // draw canvas for chart
        const chartSVG = d3.select("#chart-container")
            .append("svg")
                .attr("viewBox", [0, 0, (chartDimensions.width), (chartDimensions.height)].join(' '))
                .attr("width", "100%")
                .attr("height", "100%")
                .attr("id", "wrapper");

        // assign role for accessibility
        chartSVG.attr("role", "figure")
            .append("title")
            .text(chartTitle);

        // Add group for bounds
        chartBounds = chartSVG.append("g")
            .attr("id", "chart-bounds")
            .style("transform", `translate(${
                chartDimensions.margin.left
            }px, ${
                chartDimensions.margin.top
            }px)`);
    }

    function initTileChart({
        width = 500, // outer width, in pixels
        height = 500, // outer height, in pixels
        margin = 1, // default margins
        marginTop = margin, // top margin, in pixels
        marginBottom = margin, // left margin, in pixels
        marginLeft = margin, // left margin, in pixels
        marginRight = margin // right margin, in pixels
    }) {
        // set up global chart dimensions, including bounded dimensions
        tileChartDimensions = {
            width,
            height,
            margin: {
                top: marginTop,
                right: marginRight,
                bottom: marginBottom,
                left: marginLeft
            },
            boundedWidth: width - marginLeft - marginRight,
            boundedHeight: height - marginTop - marginBottom
        }

        tileChartBounds = chartBounds.append("g")
            .attr("id", "tile-chart-bounds")
            .style("transform", `translate(${
                tileChartDimensions.margin.left
            }px, ${
                tileChartDimensions.margin.top
            }px)`);


        // Initialize scales
        // initXScale()
        initYScale()
        initTileColorScale()

        // Initialize axes
        initYAxis({bounds: tileChartBounds})

        // Add groups for visual elements
        tileChartBounds.append("g")
            .attr("class", "rects");
        tileChartBounds.append("g")
            .attr("class", "annotations");
    }

    function initBarChart({
        width = 500, // outer width, in pixels
        height = 500, // outer height, in pixels
        margin = 1, // default margins
        marginTop = margin, // top margin, in pixels
        marginBottom = margin, // left margin, in pixels
        marginLeft = margin, // left margin, in pixels
        marginRight = margin, // right margin, in pixels
        translateX = translateX // amount to translate in x direction
    }) {
        // set up global chart dimensions, including bounded dimensions
        barChartDimensions = {
            width,
            height,
            margin: {
                top: marginTop,
                right: marginRight,
                bottom: marginBottom,
                left: marginLeft
            },
            boundedWidth: width - marginLeft - marginRight,
            boundedHeight: height - marginTop - marginBottom
        }

        barChartBounds = chartBounds.append("g")
            .attr("id", "bar-chart-bounds")
            .style("transform", `translate(${
                translateX + barChartDimensions.margin.left
            }px, ${
                barChartDimensions.margin.top
            }px)`);


        // Initialize scales
        initBarXScale()
        initBarYScale()

        // Initialize axes
        initXAxis({axis: barXAxis, bounds: barChartBounds, chartDims: barChartDimensions})

        // Add groups for visual elements
        barChartBounds.append("g")
            .attr("class", "bars");
    }

    function initScatterChart({
        width = 500, // outer width, in pixels
        height = 500, // outer height, in pixels
        margin = 1, // default margins
        marginTop = margin, // top margin, in pixels
        marginBottom = margin, // left margin, in pixels
        marginLeft = margin, // left margin, in pixels
        marginRight = margin, // right margin, in pixels
        translateX = translateX // amount to translate in x direction
    }) {
        // set up global chart dimensions, including bounded dimensions
        scatterChartDimensions = {
            width,
            height,
            margin: {
                top: marginTop,
                right: marginRight,
                bottom: marginBottom,
                left: marginLeft
            },
            boundedWidth: width - marginLeft - marginRight,
            boundedHeight: height - marginTop - marginBottom
        }

        scatterChartBounds = chartBounds.append("g")
            .attr("id", "scatter-chart-bounds")
            .style("transform", `translate(${
                translateX + scatterChartDimensions.margin.left
            }px, ${
                scatterChartDimensions.margin.top
            }px)`);


        // Initialize scales
        initScatterXScale()

        // Add groups for visual elements
        scatterChartBounds.append("g")
            .attr("class", "points");
    }

    function initBarXScale() {
        // scale for x axis (domain set in `drawBarChart()`)
        barXScale = d3.scaleLinear()
            .range([0, barChartDimensions.boundedWidth]);
    }

    function initBarYScale() {
        // scale for the y axis (domain set in `drawBarChart()`)
        barYScale = d3.scaleBand()
            .range([0, tileChartDimensions.boundedHeight]);
    }

    function initScatterXScale() {
        // scale for the x axis (domain set in `drawScatterChart()`)
        scatterXScale = d3.scaleBand()
            .range([0, scatterChartDimensions.boundedWidth]);
    }

    function initXAxis({
        axis,
        bounds,
        chartDims
    }) {
        // add group for x axis
        axis = bounds.append("g")
            .attr("id", "x-axis")
            .attr("class", "axis")
            .attr("transform", `translate(0,${chartDims.boundedHeight})`)
            .attr("aria-hidden", true); // hide from screen reader
    }

    function initYScale() {
        // scale for the y axis (domain set in `drawTileChart()`)
        yScale = d3.scaleLinear()
            .range([0, tileChartDimensions.boundedHeight]);
    }

    function initYAxis({
        bounds
    }) {
        // add group for y axis
        yAxis = bounds.append("g")
            .attr("id", "y-axis")
            .attr("class", "axis")
            .attr("aria-hidden", true);
    }

    function drawAxis({
        axis,
        axisScale,
        axisFxn,
        chartDims
    }, {
        axisTitle = '',
        titleX = -chartDims.boundedHeight / 2,
        titleY = -chartDims.margin.left,
        titleTextAnchor = "middle",
        titleBaseline = "text-before-edge",
        titleAngle = -90,
        tickSize = 0,
        tickPadding = 5,
        tickType = "numeric",
        tickFormat = ".5f",
        customSuffix = null,
        textAngle = 0,
        keepDomain = true,
    }) {
        // generate axis
        // if numeric ticks, include specification of format
        if (tickType == "numeric" && !customSuffix) {
            axis
                .call(d3[axisFxn](axisScale).tickSize(tickSize).tickPadding(tickPadding).tickFormat(d3.format(tickFormat)));
        } else if (tickType == "numeric" && customSuffix) {
            axis
                .call(d3[axisFxn](axisScale).tickSize(tickSize).tickPadding(tickPadding).tickFormat(d => d3.format(tickFormat)(d) + ' ' + customSuffix));
        } else {
            axis
                .call(d3[axisFxn](axisScale).tickSize(tickSize).tickPadding(tickPadding));
        }

        // if specified, remove axis line
        if (!keepDomain) {
            axis.select(".domain").remove();
        }

        // add class to axis labels and rotate as specified
        axis.selectAll('text')
            .attr("class", "axis-text")
            .attr("transform", `rotate(${textAngle})`);

        // add axis title
        axis
            .append("text")
            .attr("class", "axis-title")
            .attr("x", titleX)
            .attr("y", titleY)
            .attr("transform", `rotate(${titleAngle})`)
            .attr("text-anchor", titleTextAnchor)
            .attr("dominant-baseline", titleBaseline)
            .attr("role", "presentation")
            .attr("aria-hidden", true)
            .text(axisTitle);
    }

    // function drawXAxis({
    //     axisTitle = '',
    //     titleX = tileChartDimensions.boundedWidth / 2,
    //     titleY = tileChartDimensions.margin.bottom,
    //     titleTextAnchor = "middle",
    //     titleBaseline = "text-after-edge",
    //     titleAngle = 0,
    //     tickSize = 0,
    //     tickPadding = 5,
    //     tickType = 'numeric',
    //     tickFormat = ".2f",
    //     customSuffix = null,
    //     textAngle = 0, 
    //     keepDomain = true,
    // }) {
    //     drawAxis({
    //         axis: xAxis,
    //         axisScale: xScale,
    //         axisFxn: 'axisBottom'
    //     }, {
    //         axisTitle: axisTitle,
    //         titleX: titleX,
    //         titleY: titleY,
    //         titleTextAnchor: titleTextAnchor,
    //         titleBaseline: titleBaseline,
    //         titleAngle: titleAngle,tickSize: tickSize,
    //         tickPadding: tickPadding,
    //         tickType: tickType,
    //         tickFormat: tickFormat,
    //         customSuffix: customSuffix,
    //         textAngle: textAngle,
    //         keepDomain: keepDomain,
    //     })
    // }

    function drawYAxis({
        axis,
        axisScale,
        chartDims
    }, {
        axisTitle = '',
        titleX = -chartDims.boundedHeight / 2,
        titleY = -chartDims.margin.left,
        titleTextAnchor = "middle",
        titleBaseline = "text-before-edge",
        titleAngle = -90,
        tickSize = 0,
        tickPadding = 5,
        tickType = 'numeric',
        tickFormat = ".2f",
        customSuffix = null,
        textAngle = 0,
        keepDomain = true,
    }) {
        drawAxis({
            axis: axis,
            axisScale: axisScale,
            axisFxn: 'axisLeft',
            chartDims: chartDims
        }, {
            axisTitle: axisTitle,
            titleX: titleX,
            titleY: titleY,
            titleTextAnchor: titleTextAnchor,
            titleBaseline: titleBaseline,
            titleAngle: titleAngle,
            tickSize: tickSize,
            tickPadding: tickPadding,
            tickType: tickType,
            textAngle: textAngle,
            tickFormat: tickFormat,
            customSuffix: customSuffix,
            keepDomain: keepDomain,
        })
    }

    function initTileColorScale() {
        tileColorScale = d3.scaleSequential()           
            .interpolator(d3.interpolateGreys);
    }

    function initBarColorScale(data) {
        barColorScale = d3.scaleOrdinal()
            .domain(data)
            .range(data.map(item => barColors[item]));
    }

    function initScatterColorScale(data) {
        scatterColorScale = d3.scaleOrdinal()
            .domain(data)
            .range(data.map(item => scatterColors[item]));
    }

    function drawTileChart(data) {
        //////////////////////////////
        /////    PROCESS DATA    /////
        //////////////////////////////

        ///////////////////////////////////////////
        /////    SET UP ACCESSOR FUNCTIONS    /////
        ///////////////////////////////////////////
        const yAccessor = d => d.depth_cm;
        // const xAccessor = d => d.total;
        const colorAccessor = d => d.total;
        const identifierAccessor = d => 'depth-' + d.depth_cm;

        ///////////////////////////////////////////
        /////    FINISH SETTING UP Y SCALE    /////
        ///////////////////////////////////////////
        // set domain for yScale, based on data
        yScale
            .domain([0, d3.max(data, yAccessor) + 10]);
        drawYAxis(
            {axis: yAxis, axisScale: yScale, chartDims: tileChartDimensions}, 
            {tickFormat: ".0f", customSuffix: 'cm', tickSize: 3, keepDomain: false}
        )

        ///////////////////////////////////
        /////    SET UP COLOR SCALE   /////
        ///////////////////////////////////
        tileColorScale
            .domain(d3.extent(data, colorAccessor));

        ////////////////////////////////////
        /////    ADD CHART ELEMENTS    /////
        ////////////////////////////////////
        const annotationGap = tileChartDimensions.boundedWidth * 0.5;
        const annotationBuffer = annotationGap * 0.2;
        // draw chart
        tileChartBounds.select('.rects') // selects our group we set up to hold chart elements
            .selectAll(".rect") // empty selection
                .data(data) // bind data
                .enter() // instantiate chart element for each element of data
                .append("rect") // append a rectangle for each element
                    .attr("class", "rect")
                    .attr("id", d => 'bar-' + identifierAccessor(d))
                    .attr("x", annotationGap)
                    .attr("y", d => yScale(yAccessor(d)))
                    .attr("height", tileChartDimensions.boundedHeight / data.length)
                    .attr("width", tileChartDimensions.boundedWidth - annotationGap)
                    .style("fill", d => tileColorScale(colorAccessor(d)));
        // draw year bands
        tileChartBounds.select(".annotations")
            .append("rect")
                .attr("class", "year-bands")
                .attr("x", annotationGap / 2 + annotationBuffer)
                .attr("y", 0)
                .attr("height", yScale(372))
                .attr("width", 3)

        tileChartBounds.select(".annotations")
            .append("rect")
                .attr("class", "year-bands")
                .attr("x", annotationGap / 2 + annotationBuffer)
                .attr("y", yScale(378))
                .attr("height", tileChartDimensions.boundedHeight - yScale(378))
                .attr("width", 3)

        tileChartBounds.select(".annotations")
            .append("text")
                .attr("x", - yScale(372) / 2)
                .attr("y", annotationGap / 2)
                .attr("transform", "rotate(-90)")
                .attr("text-anchor", "middle")
                .text("2016 accumulation")
        
        tileChartBounds.select(".annotations")
            .append("text")
                .attr("x", - yScale(378) - ((tileChartDimensions.boundedHeight - yScale(378)) / 2))
                .attr("y", annotationGap / 2)
                .attr("transform", "rotate(-90)")
                .attr("text-anchor", "middle")
                .text("2015 accumulation")
    }

    function drawBarChart(data) {
        ///////////////////////////////////////////
        /////    SET UP ACCESSOR FUNCTIONS    /////
        ///////////////////////////////////////////
        const yAccessor = d => d.depth_cm;
        const xAccessor = d => d.picogram_per_mL;
        const keyAccessor = d => d.sugar;
        const colorAccessor = d => d.sugar;

        //////////////////////////////
        /////    PROCESS DATA    /////
        //////////////////////////////
        // sort data by bar order, so species with shared color plot together
        data.sort((a,b) => (a.bar_order > b.bar_order) ? 1 : ((b.bar_order > a.bar_order) ? -1 : 0))
        const series = d3.stack()
            .keys(d3.union(data.map(d => keyAccessor(d)))) // distinct series keys, in input order
            .value(([, D], key) => xAccessor(D.get(key))) // get value for each series key and stack
            .order(d3.stackOrderNone)
            (d3.index(data, d => yAccessor(d), d => keyAccessor(d)));

        ///////////////////////////////////////////
        /////    FINISH SETTING UP Y SCALE    /////
        ///////////////////////////////////////////
        // set domain for yScale, based on data
        barYScale
            .domain([... new Set(data.map(d => yAccessor(d)))]);
        
        ///////////////////////////////////////////
        /////    FINISH SETTING UP X SCALE    /////
        ///////////////////////////////////////////
        // // set domain for xScale
        barXScale
            .domain([0, d3.max(series, d => d3.max(d, d => d[1]))]);
        // drawXAxis({axisTitle: 'Climate vulnerability'})
        
        ///////////////////////////////////
        /////    SET UP COLOR SCALE   /////
        ///////////////////////////////////
        const colorCategories = [... new Set(data.map(colorAccessor))];
        initBarColorScale(colorCategories)

        ////////////////////////////////////
        /////    ADD CHART ELEMENTS    /////
        ////////////////////////////////////
        // draw chart
        const barGroups = barChartBounds.selectAll('.bars')
            .selectAll(".bar") // empty selection
            .data(series)
            .enter()
            .append("g")
                .attr("class", d => d.key)
        
        barGroups.selectAll(".bar") //empty selection
            .data(D => D.map(d => (d.key = D.key, d)))
            .enter()
            .append("rect")
                .attr("class", d => d.key + ' depth' + d.data[0])
                .attr("y", d => barYScale(d.data[0]))
                .attr("x", d => barXScale(d[0]))
                .attr("height", barYScale.bandwidth())
                .attr("width", d => barXScale(d[1]) - barXScale(d[0]))
                .style("fill", d => barColorScale(d.key))
    }

    function drawScatterChart(data) {
        //////////////////////////////
        /////    PROCESS DATA    /////
        //////////////////////////////

        ///////////////////////////////////////////
        /////    SET UP ACCESSOR FUNCTIONS    /////
        ///////////////////////////////////////////
        const yAccessor = d => d.depth_cm;
        const xAccessor = d => d.vegetation_type;
        const colorAccessor = d => d.vegetation_type;
        const identifierAccessor = d => 'depth-' + d.depth_cm + '-' + d.vegetation_type;

        ///////////////////////////////////////////
        /////    FINISH SETTING UP Y SCALE    /////
        ///////////////////////////////////////////
        // use bar chart y scale
        
        ///////////////////////////////////////////
        /////    FINISH SETTING UP X SCALE    /////
        ///////////////////////////////////////////
        // set domain for xScale
        scatterXScale
            .domain([... new Set(data.map(d => xAccessor(d)))]);

        ///////////////////////////////////
        /////    SET UP COLOR SCALE   /////
        ///////////////////////////////////
        const colorCategories = [... new Set(data.map(colorAccessor))];
        initScatterColorScale(colorCategories)

        ////////////////////////////////////
        /////    ADD CHART ELEMENTS    /////
        ////////////////////////////////////
        // draw chart
        scatterChartBounds.select('.points') // selects our group we set up to hold chart elements
            .selectAll(".point") // empty selection
                .data(data) // bind data
                .enter() // instantiate chart element for each element of data
                .append("circle") // append a rectangle for each element
                    .attr("class", "point")
                    .attr("id", d => 'point-' + identifierAccessor(d))
                    .attr("cx", d => scatterXScale(xAccessor(d)))
                    .attr("cy", d => barYScale(yAccessor(d)) + barYScale.bandwidth()/2)
                    .attr("r", 4)
                    .style("fill", d => scatterColorScale(colorAccessor(d)));
    }
</script>

<style scoped lang="scss">
    #wildfire-aerosols-grid-container {
        display: grid;
        max-width: 1200px;
        grid-template-columns: 10% calc(80% - 4rem) 10%;
        grid-template-rows: auto max-content;
        grid-template-areas:
            "text text text"
            "prev chart next";
        margin: 2rem auto 0 auto;
        column-gap: 2rem;
        row-gap: 3rem;
        @media only screen and (max-width: 600px) {
            width: 90vw;
            grid-template-rows: auto max-content;
            grid-template-areas:
                "chart chart chart"
                "prev text next";
        }
    }
    #chart-container {
        grid-area: chart;
        width: 100%;
    }
    #aerosol-text-container {
        grid-area: text;
        height: 10vh;
        @media screen and (max-height: 770px) {
            height: 30vh;
        }
        @media only screen and (max-width: 600px) {
            height: auto;
        }
    }
    .flip-button {
        height: 5rem;
        width: 5rem;
        align-self: center;
        border-radius: 5rem;
        border: solid 0.75px var(--medium-grey);
        cursor: pointer;
        box-shadow: 0px 0px 4px rgba(39,44,49,.3);
        @media only screen and (max-width: 600px) {
            height: 3rem;
            width: 3rem;
            align-self: start;
        }
    }
    #aerosol-prev {
        grid-area: prev;
        justify-self: end;
    }
    #aerosol-next {
        grid-area: next;
        justify-self: start;
    }
    button:hover:after {
        top: 0px;
        left: 0px;
    }
    button:hover {
        box-shadow: rgba(39,44,49,.7) 2px 2px 4px -2px;
        transform: translate3d(0, 2px, 0);
    }
    button:disabled {
        background-color: #b4b2b2;
        cursor: not-allowed;
        border-color: #b4b2b2;
        color: #b4b2b2;
    }
    button:disabled:hover {
        box-shadow: 0px 0px 4px rgba(39,44,49,.3);
        transform: translate3d(0, 0px, 0);
    }
    button:disabled:after {
        background-color: #b4b2b2;
    }
    .fa {
        color: var(--color-text);
        opacity: 1;
        @media only screen and (max-width: 600px) {
            font-size: 1.6rem;
        }
    }
    button:disabled .fa {
        opacity: 0.2;
    }
</style>
<style lang="scss">
/* css for elements added/classed w/ d3 */
    .axis-text {
        font-size: 1.6rem;
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
    .year-bands {
        fill: var(--medium-grey);
    }
</style>