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
                <button id="aerosol-prev" class="flip-button" @click="currentIndex--; clicked()" :disabled="isFirstImage || justClicked">
                    <font-awesome-icon :icon="{ prefix: 'fas', iconName: 'arrow-left' }"  class="fa fa-arrow-left"/>
                </button>
                <button id="aerosol-next" class="flip-button" @click="currentIndex++; clicked()" :disabled="isLastImage || justClicked">
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
    import { computed, onMounted, ref, watch } from "vue";
    import { isMobile } from 'mobile-device-detect';
    import * as d3 from 'd3';
    import VizSection from '@/components/VizSection.vue';

    // define props
    const props = defineProps({
        text: { type: Object }
    })

    // global variables
    const publicPath = import.meta.env.BASE_URL;
    const mobileView = isMobile;
    const tileDataFile = 'fii_core4particulates.csv';
    const barDataFile = 'fii_core4sugars.csv';
    const scatterDataFile = 'fii_core4biomass.csv';
    const tileData = ref();
    const barData = ref();
    const scatterData = ref();
    const currentIndex = ref(1);
    const justClicked = ref(false);
    const nIndices = 3;
    const chart = ref(null);
    let chartSVG;
    const chartTitle = 'Title of chart';
    const chartHeight = mobileView ? window.innerHeight * 0.6 : window.innerHeight * 0.7;
    let chartWidth;
    let chartDimensions;
    let chartBounds;
    let tileChartTranslateX1;
    let tileChartTranslateX2;
    let tileChartTranslateX3;
    let tileChartDimensions;
    let tileChartBounds;
    let tileColorScale;
    let yScale;
    let yAxis;
    let barChartTranslateX2;
    let barChartTranslateX3;
    let barChartDimensions;
    let barChartBounds;
    let barXScale;
    let barYScale;
    let barXAxis;
    const barXAxisPosition = 'top';
    let barColorCategories;
    const barColors = {Mannosan: '#000000', Galactosan: '#989898', Levoglucosan: '#c8c8c8'};
    let barColorScale;
    let scatterChartTranslateX3;
    let scatterChartDimensions;
    let scatterChartBounds;
    let scatterXScale;
    let scatterColorCategories;
    const scatterColors = {grass: '#c49051', hardwood: '#3c475a', softwood: '#729C9D'};
    let scatterColorScale;
    const transitionLength = 20;

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

    watch(currentIndex, () => {
        updateChartView(currentIndex.value)
    })

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
                    width: chartWidth,
                    height: chartHeight,
                    margin: 0
                })

                const defaultMargin = mobileView ? 5 : 10;
                const sharedTopMargin = mobileView ? 100 : 120;
                const sharedBottomMargin = mobileView ? 0 : 10;

                const tileChartWidth = chartDimensions.boundedWidth / 3
                tileChartTranslateX1 = mobileView ? 40 : 350;
                tileChartTranslateX2 = mobileView ? 40 : 200;
                tileChartTranslateX3 = mobileView ? 75 : 80;
                initTileChart({
                    width: tileChartWidth,
                    height: chartHeight,
                    margin: defaultMargin,
                    marginLeft: mobileView ? 55: 80,
                    marginRight: mobileView ? 5: 40,
                    marginTop: sharedTopMargin,
                    marginBottom: sharedBottomMargin,
                    translateX: tileChartTranslateX1
                });

                const barChartWidth = chartDimensions.boundedWidth / 3
                barChartTranslateX2 = mobileView ? tileChartWidth + 40 : tileChartWidth + 150;
                barChartTranslateX3 = mobileView ? tileChartWidth + 10 : tileChartWidth + 80;
                initBarChart({
                    width: barChartWidth,
                    height: chartHeight,
                    margin: defaultMargin,
                    marginLeft: mobileView ? 20 : 80,
                    marginTop: sharedTopMargin,
                    marginBottom: sharedBottomMargin,
                    translateX: barChartTranslateX2});

                const scatterChartWidth = chartDimensions.boundedWidth / 3
                scatterChartTranslateX3 = tileChartWidth + barChartWidth + 80;
                initScatterChart({
                    width: scatterChartWidth,
                    height: chartHeight,
                    margin: defaultMargin,
                    marginLeft: mobileView ? 20 : scatterChartWidth / 3,
                    marginRight: scatterChartWidth / 3,
                    marginTop: sharedTopMargin,
                    marginBottom: sharedBottomMargin,
                    translateX: scatterChartTranslateX3});

                // draw charts, hiding bar and scatter charts to start
                drawTileChart(tileData.value);
                drawBarChart(barData.value);
                barChartBounds
                    .attr("display", "none");
                drawScatterChart(scatterData.value);
                scatterChartBounds
                    .attr("display", "none");

                // add legends
                addTileLegend()
                addBarLegend()
                addScatterLegend()
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

    function clicked() {
        justClicked.value = true;
        setTimeout(function(){justClicked.value = false;}, transitionLength);
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
        chartSVG = d3.select("#chart-container")
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
        marginRight = margin, // right margin, in pixels
        translateX // starting x translation
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
                translateX
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
        initBarXAxis(
            {
                bounds: barChartBounds, 
                chartDims: barChartDimensions,
                axisPosition: barXAxisPosition
            }
        )

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

    function initBarXAxis({
        bounds,
        chartDims,
        axisPosition = 'bottom'
    }) {
        // add group for x axis
        barXAxis = bounds.append("g")
            .attr("id", "x-axis")
            .attr("class", "axis")
            .attr("transform", `translate(0,${axisPosition === 'bottom' ? chartDims.boundedHeight : 0})`)
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
        axisSubtitle = null,
        titleX = -chartDims.boundedHeight / 2,
        titleY = -chartDims.margin.left,
        titleTextAnchor = "middle",
        titleBaseline = "text-before-edge",
        titleAngle = -90,
        nticks = null,
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
                .call(d3[axisFxn](axisScale).ticks(nticks).tickSize(tickSize).tickPadding(tickPadding).tickFormat(d3.format(tickFormat)));
        } else if (tickType == "numeric" && customSuffix) {
            axis
                .call(d3[axisFxn](axisScale).ticks(nticks).tickSize(tickSize).tickPadding(tickPadding).tickFormat(d => d3.format(tickFormat)(d) + ' ' + customSuffix));
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
        axisTitle = axis
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

        if (axisSubtitle) {
            axisTitle.append("tspan")
                .attr("class", "axis-subtitle")
                .attr("x", titleX)
                .attr("y", titleY)
                .attr("dy", "1em")
                .attr("transform", `rotate(${titleAngle})`)
                .attr("text-anchor", titleTextAnchor)
                .attr("dominant-baseline", titleBaseline)
                .attr("role", "presentation")
                .attr("aria-hidden", true)
                .text(axisSubtitle);
        }
    }

    function drawXAxis({
        axis,
        axisScale,
        chartDims
    }, {
        axisPosition = 'bottom',
        axisTitle = '',
        axisSubtitle = null,
        titleX = chartDims.boundedWidth / 2,
        titleY = axisPosition === 'bottom' ? chartDims.margin.bottom : -chartDims.margin.top,
        titleTextAnchor = "middle",
        titleBaseline = axisPosition === 'bottom' ? "text-after-edge" : "text-before-edge",
        titleAngle = 0,
        nticks = null,
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
            axisFxn: axisPosition === 'bottom' ? 'axisBottom' : 'axisTop'
        }, {
            axisTitle: axisTitle,
            axisSubtitle: axisSubtitle,
            titleX: titleX,
            titleY: titleY,
            titleTextAnchor: titleTextAnchor,
            titleBaseline: titleBaseline,
            titleAngle: titleAngle,
            nticks: nticks,
            tickSize: tickSize,
            tickPadding: tickPadding,
            tickType: tickType,
            tickFormat: tickFormat,
            customSuffix: customSuffix,
            textAngle: textAngle,
            keepDomain: keepDomain,
        })
    }

    function drawYAxis({
        axis,
        axisScale,
        chartDims
    }, {
        axisTitle = '',
        axisSubtitle = null,
        titleX = -chartDims.boundedHeight / 2,
        titleY = -chartDims.margin.left,
        titleTextAnchor = "middle",
        titleBaseline = "text-before-edge",
        titleAngle = -90,
        nticks = null,
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
            axisSubtitle: axisSubtitle,
            titleX: titleX,
            titleY: titleY,
            titleTextAnchor: titleTextAnchor,
            titleBaseline: titleBaseline,
            titleAngle: titleAngle,
            nticks: nticks,
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
            {
                axis: yAxis, 
                axisScale: yScale, 
                chartDims: tileChartDimensions
            }, 
            {
                tickFormat: ".0f", 
                customSuffix: 'cm', 
                tickSize: 3, 
                keepDomain: false
            }
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

    function addTileLegend() {
        // build list of posible counts (0 to 366)
        let count_list = [];
        for (let i = 1; i <= tileColorScale.domain()[1]; i++) {
            count_list.push(i);
        }

        // define gradient for legend
        chartSVG.append("defs")
            .append("linearGradient")
            .attr("id", "gradient-particles")
            .attr("x1", "0%").attr("y1", "0%")
            .attr("x2", "100%").attr("y2", "0%")
            .selectAll("stop")
            .data(count_list)
                .enter()
                .append("stop")
                .attr("offset", (d,i) => i/(count_list.length-1))
                .attr("stop-color", d => tileColorScale(d))

        const legendGroup = tileChartBounds.append("g")
            .attr("id", "tile-chart-legend")

        // append legend title
        legendGroup.append("text")
              .attr("class", "axis-title")
              .attr("text-anchor", "start")
              .attr("x", 0)
              .attr("y", -tileChartDimensions.margin.top)
              .attr("dominant-baseline", "text-before-edge")
              .text("Particulate count")

        // append legend text
        
        const xBuffer = 5;
        legendGroup.append("text")
              .attr("class", "axis-subtitle")
              .attr("text-anchor", "end")
              .attr("dominant-baseline", "central")
              .attr("x", - xBuffer)
              .attr("y", -tileChartDimensions.margin.top / 2)
              .text("low")

        legendGroup.append("text")
              .attr("class", "axis-subtitle")
              .attr("text-anchor", "start")
              .attr("dominant-baseline", "central")
              .attr("x", tileChartDimensions.boundedWidth / 2 + xBuffer)
              .attr("y", -tileChartDimensions.margin.top / 2)
              .text("high")
        
        // append legend rectangle
        const rectHeight = tileChartDimensions.margin.top / 4
        legendGroup.append("rect")
              .attr("class", "c1p2 matrixLegend")
              .attr("width", tileChartDimensions.boundedWidth / 2)
              .attr("height", rectHeight)
              .attr("fill", "url(#gradient-particles)")
              .attr("x", 0)
              .attr("y", -tileChartDimensions.margin.top / 2 - rectHeight / 2)
    }

    function drawBarChart(data) {
        ///////////////////////////////////////////
        /////    SET UP ACCESSOR FUNCTIONS    /////
        ///////////////////////////////////////////
        const yAccessor = d => d.depth_cm;
        const xAccessor = d => d.picogram_per_mL/1000; // convert to nanogram
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
        drawXAxis(
            {
                axis: barXAxis, 
                axisScale: barXScale, 
                chartDims: barChartDimensions
            }, 
            {
                axisPosition: barXAxisPosition, 
                axisTitle: 'Sugar concentration', 
                axisSubtitle: 'nanogram/mL',
                nticks: 4,
                tickFormat: '.0f', 
                tickSize: 3
            }
        )
        
        ///////////////////////////////////
        /////    SET UP COLOR SCALE   /////
        ///////////////////////////////////
        barColorCategories = [... new Set(data.map(colorAccessor))];
        initBarColorScale(barColorCategories)

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

    function addBarLegend() {
        const barLegendGroup = barChartBounds.append("g")
            .attr("id", "bar-chart-legend")

        // // Add legend title
        // barLegendGroup.append("text")
        //     .text('Sugars')
        //     .attr("id", "legend-title")
        //     .attr("class", "axis-title")
        //     .attr("y", barChartDimensions.margin.top / 2)
        //     .attr("dominant-baseline", "central")

        const legendRectSize = barYScale.bandwidth();
        const interItemSpacing = mobileView ? 15 : 10;
        const intraItemSpacing = 6;

        // Append group for each legend entry
        const legendGroup = barLegendGroup.selectAll(".legend-item")
            .data(barColorCategories)
            .enter()
            .append("g")
            .attr("class", "legend-item")

        // Add rectangles for each group
        legendGroup.append("rect")
            .attr("class", "legend-rect")
            .attr("x", 0)
            .attr("y", -barChartDimensions.margin.top / 2 - legendRectSize / 2.5)
            .attr("width", legendRectSize)
            .attr("height", legendRectSize)
            .style("fill", d => barColorScale(d))
        
        // Add text for each group
        legendGroup.append("text")
            .attr("class", "legend-text")
            .attr("x", legendRectSize + intraItemSpacing) // put text to the right of the rectangle
            .attr("y", -barChartDimensions.margin.top / 2)
            .attr("text-anchor", "start") // left-align text
            .attr("dominant-baseline", "central")
            .text(d => d);

        // Position legend groups
        // https://stackoverflow.com/questions/20224611/d3-position-text-element-dependent-on-length-of-element-before
        const xBuffer = 6; // set xBuffer for use in mobile row x translations
        legendGroup
            .attr("transform", (d, i) => {
            // Compute total width of preceeding legend items, with spacing
            // Start with width of legend title
            const titleWidth = 0 //d3.select('#legend-title')._groups[0][0].getBBox().width + interItemSpacing;
            
            // Begin right of legend
            let cumulativeWidth = titleWidth;
            // if (mobileView) {
                // On mobile, only use preceding items in same row to find cumulative width
                // row 1: items 0 and 1
                if (i < 2) {
                for (let j = 0; j < i; j++) {
                    cumulativeWidth = cumulativeWidth + legendGroup._groups[0][j].getBBox().width + interItemSpacing;
                }
                }
                // row 2: items 2, 3 and 4
                else if (i < 5) {
                for (let j = 2; j < i; j++) {
                    cumulativeWidth = cumulativeWidth + legendGroup._groups[0][j].getBBox().width + interItemSpacing;
                }
                }
                // row 3: item 5 
                else if (i === 5) {
                for (let j = 2; j < i; j++) {
                    // Actually storing width of row 2 here, to use to set selfSupplyEnd
                    cumulativeWidth = cumulativeWidth + legendGroup._groups[0][j].getBBox().width + interItemSpacing;
                }
                }
            // } else {
                // on desktop, iterate through all preceding items to find cumulative width, since all items in 1 row
                // for (let j = 0; j < i; j++) {
                // cumulativeWidth = cumulativeWidth + legendGroup._groups[0][j].getBBox().width + interItemSpacing;
                // }
            // }

            let yTranslation = 0;
            // Determine x and y translation
            // set y translation for each row
            // adjust row starting position for 2nd and third rows by -titleWidth
            // if (mobileView) {
                if (i < 2) {
                    yTranslation = 0;
                } else if (i < 5) {
                    yTranslation = mobileView ? legendRectSize * 4 : legendRectSize * 2;
                    cumulativeWidth = cumulativeWidth - titleWidth;
                } else {
                    yTranslation = legendRectSize * 5.75
                    cumulativeWidth = xBuffer; // for last item just translate by xBuffer
                } 
            // }

            // translate each group by that width and height
            return "translate(" + cumulativeWidth + "," + yTranslation + ")"
        })
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
        scatterColorCategories = [... new Set(data.map(colorAccessor))];
        initScatterColorScale(scatterColorCategories)

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

    function addScatterLegend() {
        const scatterLegendGroup = scatterChartBounds.append("g")
            .attr("id", "bar-chart-legend")

        // Add legend title
        scatterLegendGroup.append("text")
            .attr("id", "legend-title")
            .attr("class", "axis-title")
            .attr("x", scatterChartDimensions.boundedWidth / 2)
            .attr("y", -scatterChartDimensions.margin.top)
            .attr("text-anchor", "middle")
            .attr("dominant-baseline", "text-before-edge")
            .attr("text-width", scatterChartDimensions.boundedWidth + (scatterChartDimensions.margin.left + scatterChartDimensions.margin.right) / 2)
            .text('Burned vegetation type')
            .call(d => wrap(d))

        const legendPointSize = 4;
        const interItemSpacing = mobileView ? 15 : 10;
        const intraItemSpacing = 6;

        // Append group for each legend entry
        const legendGroup = scatterLegendGroup
            .append("g")
            .style("transform", `translate(${
                -scatterChartDimensions.margin.left / 4
            }px, 0px)`);

        const legendGroups = legendGroup.selectAll(".legend-item")
            .data(scatterColorCategories)
            .enter()
            .append("g")
            .attr("class", "legend-item")

        // Add rectangles for each group
        legendGroups.append("circle")
            .attr("class", "legend-point")
            .attr("cx", 0)
            .attr("cy", -scatterChartDimensions.margin.top / 2 - legendPointSize / 2.5)
            .attr("r", legendPointSize)
            .style("fill", d => scatterColorScale(d))
        
        // Add text for each group
        legendGroups.append("text")
            .attr("class", "legend-text")
            .attr("x", legendPointSize + intraItemSpacing) // put text to the right of the rectangle
            .attr("y", -scatterChartDimensions.margin.top / 2)
            .attr("text-anchor", "start") // left-align text
            .attr("dominant-baseline", "central")
            .text(d => d);

        // Position legend groups
        // https://stackoverflow.com/questions/20224611/d3-position-text-element-dependent-on-length-of-element-before
        const xBuffer = 6; // set xBuffer for use in mobile row x translations
        legendGroups
            .attr("transform", (d, i) => {
                // Compute total width of preceeding legend items, with spacing
                // Start with width of legend title
                const titleWidth = 0 //d3.select('#legend-title')._groups[0][0].getBBox().width + interItemSpacing;
                
                // Begin right of legend
                let cumulativeWidth = titleWidth;
                // if (mobileView) {
                    // On mobile, only use preceding items in same row to find cumulative width
                    // row 1: items 0 and 1
                    if (i < 2) {
                    for (let j = 0; j < i; j++) {
                        cumulativeWidth = cumulativeWidth + legendGroups._groups[0][j].getBBox().width + interItemSpacing;
                    }
                    }
                    // row 2: items 2, 3 and 4
                    else if (i < 5) {
                    for (let j = 2; j < i; j++) {
                        cumulativeWidth = cumulativeWidth + legendGroups._groups[0][j].getBBox().width + interItemSpacing;
                    }
                    }
                    // row 3: item 5 
                    else if (i === 5) {
                    for (let j = 2; j < i; j++) {
                        // Actually storing width of row 2 here, to use to set selfSupplyEnd
                        cumulativeWidth = cumulativeWidth + legendGroups._groups[0][j].getBBox().width + interItemSpacing;
                    }
                    }
                // } else {
                    // on desktop, iterate through all preceding items to find cumulative width, since all items in 1 row
                    // for (let j = 0; j < i; j++) {
                    // cumulativeWidth = cumulativeWidth + legendGroups._groups[0][j].getBBox().width + interItemSpacing;
                    // }
                // }

                let yTranslation = 0;
                // Determine x and y translation
                // set y translation for each row
                // adjust row starting position for 2nd and third rows by -titleWidth
                // if (mobileView) {
                    if (i < 2) {
                        yTranslation = 0;
                    } else if (i < 5) {
                        yTranslation = mobileView ? barYScale.bandwidth() * 4 : barYScale.bandwidth() * 2;
                        cumulativeWidth = cumulativeWidth - titleWidth;
                    } else {
                        yTranslation = legendPointSize * 5.75
                        cumulativeWidth = xBuffer; // for last item just translate by xBuffer
                    } 
                // }

                // translate each group by that width and height
                return "translate(" + cumulativeWidth + "," + yTranslation + ")"
            })
    }

    function showChart(hiddenChartBounds, translateX, translateY) {
        hiddenChartBounds
            .attr("opacity", 0)
            .attr("display", "auto")
            .transition()
            .duration(transitionLength)
            .attr("opacity", 1)
            .style("transform", `translate(${
                translateX
            }px, ${
                translateY
            }px)`);
    }

    function hideChart(visibleChartBounds) {
        visibleChartBounds
            .transition()
            .duration(transitionLength)
            .attr("opacity", 0)
        setTimeout(function() {
            visibleChartBounds
                .attr("display", "none")
        }, transitionLength)
    }

    function moveChart(visibleChartBounds, translateX, translateY) {
        visibleChartBounds
            .transition()
            .duration(transitionLength)
            .style("transform", `translate(${
                translateX
            }px, ${
                translateY
            }px)`);
    }

    function updateChartView(index) {
        const barChartHidden = barChartBounds.node().attributes.display.value == 'none';
        const scatterChartHidden = scatterChartBounds.node().attributes.display.value == 'none';
        if (index == 1) {
            moveChart(tileChartBounds, tileChartTranslateX1, tileChartDimensions.margin.top)
            if (!barChartHidden) hideChart(barChartBounds)            
        } else if (index == 2) {
            moveChart(tileChartBounds, tileChartTranslateX2, tileChartDimensions.margin.top)
            if (barChartHidden)  {
                showChart(barChartBounds, barChartTranslateX2, barChartDimensions.margin.top)
            } else {
                moveChart(barChartBounds, barChartTranslateX2, barChartDimensions.margin.top)
            }
            if (!scatterChartHidden) hideChart(scatterChartBounds) 
        } else if (index == 3) {
            moveChart(tileChartBounds, tileChartTranslateX3, tileChartDimensions.margin.top)
            moveChart(barChartBounds, barChartTranslateX3, barChartDimensions.margin.top)
            showChart(scatterChartBounds, scatterChartTranslateX3, scatterChartDimensions.margin.top)
        }
    }

    // https://gist.github.com/mbostock/7555321
    function wrap(text) {
        text.each(function() {
            var text = d3.select(this),
            words = text.text().split(/\s|-+/).reverse(),
            word,
            line = [],
            lineNumber = 0,
            lineHeight = 1.1, // ems
            width = text.attr("text-width"),
            baseline = text.attr("dominant-baseline"),
            x = text.attr("x"),
            y = text.attr("y"),
            dy = parseFloat(text.attr("dy")),
            dx = parseFloat(text.attr("dx")),
            tspan = text.text(null).append("tspan").attr("y", y).attr("dy", dy + "em");

            while ((word = words.pop())) {
            line.push(word);
            tspan.text(line.join(" "));
                if (tspan.node().getComputedTextLength() > width) {
                line.pop();
                tspan.text(line.join(" "));
                line = [word];
                tspan = text.append("tspan").attr("x", x).attr("y", y).attr("dx", dx).attr("dy", ++lineNumber * lineHeight + dy + "em").text(word);
                }
            }

            // https://stackoverflow.com/questions/60558291/wrapping-and-vertically-centering-text-using-d3-js
            if (lineNumber > 0) {
                let lineHeightFactor;
                switch(baseline) {
                    case 'central':
                        lineHeightFactor = 0;
                        break;
                    case 'text-before-edge':
                        lineHeightFactor = 0;
                        break;
                    default:
                        lineHeightFactor = 0;
                }
                const startDy = -(lineNumber * (lineHeight / 2)) * lineHeightFactor;
                text
                    .selectAll("tspan")
                    .attr("dy", (d, i) => startDy + lineHeight * i + "em");
            }
        }
    )};

</script>

<style scoped lang="scss">
    #wildfire-aerosols-grid-container {
        display: grid;
        max-width: 1200px;
        grid-template-columns: 10% calc(80% - 4rem) 10%;
        grid-template-rows: auto max-content;
        grid-template-areas:
            "prev text next"
            "chart chart chart";
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
        align-self: start;
    }
    #aerosol-next {
        grid-area: next;
        justify-self: start;
        align-self: start;
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
        font-weight: 300;
        font-family: var(--default-font);
        user-select: none;
        @media only screen and (max-width: 600px) {
            font-size: 1.4rem;
        }
    }
    .axis-title {
        font-size: 1.8rem;
        font-family: var(--default-font);
        font-weight: 900;
        fill: var(--color-text);
        user-select: none;
        @media only screen and (max-width: 600px) {
            font-size: 1.6rem;
        }
    }
    .axis-subtitle {
        font-size: 1.8rem;
        font-family: var(--default-font);
        font-weight: 300;
        fill: var(--color-text);
        user-select: none;
        @media only screen and (max-width: 600px) {
            font-size: 1.6rem;
        }
    }
    .year-bands {
        fill: var(--medium-grey);
    }
    .legend-text {
        font-size: 1.6rem;
        font-family: var(--default-font);
        user-select: none;
        @media only screen and (max-width: 600px) {
            font-size: 1.4rem;
        }
    }
</style>