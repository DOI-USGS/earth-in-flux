<template>
    <section>
        <!---VizSection-->
        <VizSection
            id="cross-section"
            :figures="true"
            :fig-caption="false"
        >
            <template #figures>
                <div id="wildfire-aerosols-grid-container">
                    <div id="aerosol-text-container" class="text-container">
                        <p v-html="currentText" />
                    </div>
                    <button id="aerosol-next-upper" class="flip-button" @click="currentIndex++; clicked()" :disabled="isLastImage || justClicked" aria-label="Upper button to see next slide">
                        <font-awesome-icon :icon="{ prefix: 'fas', iconName: 'arrow-right' }"  class="fa fa-arrow-right"/>
                    </button>
                    <button id="aerosol-prev-upper" class="flip-button" @click="currentIndex--; clicked()" :disabled="isFirstImage || justClicked" aria-label="Upper button to see previous slide">
                        <font-awesome-icon :icon="{ prefix: 'fas', iconName: 'arrow-left' }"  class="fa fa-arrow-left"/>
                    </button>
                    <div id="chart-container" ref="chart"></div>
                    <button v-if="!mobileView" id="aerosol-next-lower" class="flip-button" @click="currentIndex++; clicked()" :disabled="isLastImage || justClicked" aria-label="Lower button to see next slide">
                        <font-awesome-icon :icon="{ prefix: 'fas', iconName: 'arrow-right' }"  class="fa fa-arrow-right"/>
                    </button>
                    <button v-if="!mobileView" id="aerosol-prev-lower" class="flip-button" @click="currentIndex--; clicked()" :disabled="isFirstImage || justClicked" aria-label="Lower button to see previous slide">
                        <font-awesome-icon :icon="{ prefix: 'fas', iconName: 'arrow-left' }"  class="fa fa-arrow-left"/>
                    </button>
                </div>
            </template>
        </VizSection>
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
            <template #aboveExplanation>
                <p v-html="text.explanation1" />
                <p v-html="text.explanation2" />
            </template>
        </VizSection>
    </section>
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
    const nIndices = 4;
    const chart = ref(null);
    let chartSVG;
    const chartTitle = 'Series of charts depicting particulate counts and the presence of wildfire biomarkers in layers of a 780-centimeter snow core';
    let chartHeight;
    let chartWidth;
    let chartDimensions;
    let chartBounds;
    let chartGap;
    let maskingRect;
    let annotationGap;
    let tileChartTranslateX1;
    let tileChartTranslateX2;
    let tileChartTranslateX3;
    let tileChartDimensions;
    let tileChartWrapper;
    let tileChartBounds;
    let tileColorScale;
    let yScale;
    let yAxis;
    let barChartTranslateX2;
    let barChartTranslateX3;
    let barChartDimensions;
    let barChartWrapper;
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
    let scatterChartWrapper;
    let scatterChartBounds;
    let scatterXScale;
    let scatterXAxis;
    const scatterXAxisPosition = 'top';
    let scatterColorCategories;
    const scatterColors = {hardwood: '#c49051', softwood: '#729C9D'};
    let scatterColorScale;
    const transitionLength = 1000;

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
                // on desktop, don't let chart height exceed 800px
                const desktopHeight = window.innerHeight < 770 ? window.innerHeight * 1.05 : Math.min(window.innerHeight * 0.75, 800);
                chartHeight = mobileView ? window.innerHeight * 0.6 : desktopHeight;
                chartWidth = chart.value.offsetWidth;
                initChart({
                    width: chartWidth,
                    height: chartHeight,
                    margin: mobileView ? 5 : 5,
                    marginLeft: mobileView ? 5 : 30
                })

                const defaultMargin = mobileView ? 5 : 10;
                const sharedTopMargin = mobileView ? 135 : 165;
                const sharedBottomMargin = mobileView ? 0 : 10;

                chartGap = mobileView ? chartDimensions.boundedWidth / 11 : chartDimensions.boundedWidth / 11;
                const tileChartWidth = mobileView ? chartGap * 3 : chartGap * 2.5;
                const barChartWidth = mobileView ? chartGap * 3 : chartGap * 3;
                const scatterChartWidth = mobileView ? chartGap * 2 : chartGap * 2;

                tileChartTranslateX1 = mobileView ? (chartDimensions.boundedWidth - tileChartWidth) / 2: (chartDimensions.boundedWidth - tileChartWidth) / 2;
                tileChartTranslateX2 = mobileView ? tileChartTranslateX1 - barChartWidth / 1.75 : tileChartTranslateX1 - barChartWidth / 1.75;
                tileChartTranslateX3 = mobileView ? tileChartTranslateX2 - scatterChartWidth + chartGap : tileChartTranslateX2 - scatterChartWidth;
                initTileChart({
                    width: tileChartWidth,
                    height: chartHeight,
                    margin: defaultMargin,
                    marginLeft: mobileView ? 5: 5,
                    marginRight: mobileView ? 5: 5,
                    marginTop: sharedTopMargin,
                    marginBottom: sharedBottomMargin,
                    translateX: tileChartTranslateX1
                });

                barChartTranslateX2 = mobileView ? tileChartTranslateX2 + tileChartWidth + chartGap : tileChartTranslateX2 + tileChartWidth + chartGap;
                barChartTranslateX3 = mobileView ? tileChartTranslateX3 + tileChartWidth + chartGap * 0.75 : tileChartTranslateX3 + tileChartWidth + chartGap;
                initBarChart({
                    width: barChartWidth,
                    height: chartHeight,
                    margin: defaultMargin,
                    marginLeft: mobileView ? 5 : 5,
                    marginRight: mobileView ? 5 : 5,
                    marginTop: sharedTopMargin,
                    marginBottom: sharedBottomMargin,
                    translateX: barChartTranslateX2});

                scatterChartTranslateX3 = mobileView ? tileChartTranslateX3 + tileChartWidth + barChartWidth + chartGap * 1.5 : tileChartTranslateX3 + tileChartWidth + barChartWidth + chartGap * 2;
                initScatterChart({
                    width: scatterChartWidth,
                    height: chartHeight,
                    margin: defaultMargin,
                    marginLeft: mobileView ? 5 : 5,
                    marginRight: mobileView ? 5 : 5,
                    marginTop: sharedTopMargin,
                    marginBottom: sharedBottomMargin,
                    translateX: scatterChartTranslateX3});

                // draw charts, hiding bar and scatter charts to start
                drawTileChart(tileData.value);
                drawBarChart(barData.value);
                barChartWrapper
                    .attr("visibility", "hidden");
                maskingRect = chartSVG
                    .append("rect")
                    .attr("id", "masking-rect")
                    .attr("x", 0)
                    .attr("y", 0)
                    .attr("width", barChartDimensions.width + scatterChartDimensions.width)
                    .attr("height", chartHeight)
                    .style("transform", `translate(${
                        tileChartTranslateX1 + tileChartDimensions.width + chartGap
                    }px, 0px)`);
                drawScatterChart(scatterData.value, 'softwood');
                scatterChartWrapper
                    .attr("visibility", "hidden");

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

        tileChartWrapper = chartBounds.append("g")
            .attr("id", "tile-chart-wrapper")
            .style("transform", `translate(${
                translateX
            }px, ${
                0
            }px)`);

        tileChartBounds = tileChartWrapper.append("g")
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
            .attr("class", "annotations");
        tileChartBounds.append("g")
            .attr("class", "rects");
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

        barChartWrapper = chartBounds.append("g")
            .attr("id", "bar-chart-wrapper")
            .style("transform", `translate(${
                translateX
            }px, ${
                0
            }px)`);

        barChartBounds = barChartWrapper.append("g")
            .attr("id", "bar-chart-bounds")
            .style("transform", `translate(${
                barChartDimensions.margin.left
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

        scatterChartWrapper = chartBounds.append("g")
            .attr("id", "scatter-chart-wrapper")
            .style("transform", `translate(${
                translateX
            }px, ${
                0
            }px)`);

        scatterChartBounds = scatterChartWrapper.append("g")
            .attr("id", "scatter-chart-bounds")
            .style("transform", `translate(${
                scatterChartDimensions.margin.left
            }px, ${
                scatterChartDimensions.margin.top
            }px)`);


        // Initialize scales
        initScatterXScale()

        // Initialize axes
        initScatterXAxis(
            {
                bounds: scatterChartBounds, 
                chartDims: scatterChartDimensions,
                axisPosition: scatterXAxisPosition
            }
        )

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

    function initScatterXAxis({
        bounds,
        chartDims,
        axisPosition = 'bottom'
    }) {
        // add group for x axis
        scatterXAxis = bounds.append("g")
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
        titleBaseline = "hanging",
        titleAngle = -90,
        titleWidth = chartDims.boundedWidth,
        wrapTitle = false,
        nticks = null,
        tickSize = 0,
        tickPadding = 5,
        tickType = "numeric",
        tickFormat = ".5f",
        customSuffix = null,
        textAngle = 0,
        keepDomain = true,
        keepLabels = true
    }) {
        // generate axis
        // if numeric ticks, include specification of format
        if (tickType == "numeric" && !customSuffix) {
            axis
                .call(d3[axisFxn](axisScale).ticks(nticks).tickSize(tickSize).tickPadding(tickPadding).tickFormat(d3.format(tickFormat)));
        } else if (tickType == "numeric" && customSuffix) {
            axis
                .call(d3[axisFxn](axisScale).ticks(nticks).tickSize(tickSize).tickPadding(tickPadding).tickFormat(d => d3.format(tickFormat)(d) + ' ' + customSuffix));
        } else if (!keepLabels) {
            axis
                .call(d3[axisFxn](axisScale).tickValues([]));
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
            .attr("dx", 0)
            .attr("dy", 0)
            .attr("transform", `rotate(${titleAngle})`)
            .attr("text-anchor", titleTextAnchor)
            .attr("dominant-baseline", titleBaseline)
            .attr("text-width", titleWidth)
            .attr("role", "presentation")
            .attr("aria-hidden", true)
            .text(axisTitle)
            .call(d => wrapTitle ? wrap(d, {shift: false}) : d);

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
        titleBaseline = axisPosition === 'bottom' ? "ideographic" : "hanging",
        titleAngle = 0,
        titleWidth = chartDims.boundedWidth,
        wrapTitle = false,
        nticks = null,
        tickSize = 0,
        tickPadding = 5,
        tickType = 'numeric',
        tickFormat = ".2f",
        customSuffix = null,
        textAngle = 0, 
        keepDomain = true,
        keepLabels = true
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
            titleWidth: titleWidth,
            wrapTitle: wrapTitle,
            nticks: nticks,
            tickSize: tickSize,
            tickPadding: tickPadding,
            tickType: tickType,
            tickFormat: tickFormat,
            customSuffix: customSuffix,
            textAngle: textAngle,
            keepDomain: keepDomain,
            keepLabels: keepLabels
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
        titleBaseline = "hanging",
        titleAngle = -90,
        titleWidth = chartDims.boundedWidth,
        wrapTitle = false,
        nticks = null,
        tickSize = 0,
        tickPadding = 5,
        tickType = 'numeric',
        tickFormat = ".2f",
        customSuffix = null,
        textAngle = 0,
        keepDomain = true,
        keepLabels = true
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
            titleWidth: titleWidth,
            wrapTitle: wrapTitle,
            nticks: nticks,
            tickSize: tickSize,
            tickPadding: tickPadding,
            tickType: tickType,
            textAngle: textAngle,
            tickFormat: tickFormat,
            customSuffix: customSuffix,
            keepDomain: keepDomain,
            keepLabels: keepLabels
        })
    }

    function initTileColorScale() {
        tileColorScale = d3.scaleSequential()           
            .range(["#efefef" ,"#000000"]);
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
                tickSize: -(chartDimensions.boundedWidth - tileChartDimensions.margin.left), 
                nticks: 20,
                keepDomain: false
            }
        )
        tileChartBounds.selectAll(".tick line")
            .attr("class", "tick-line")

        ///////////////////////////////////
        /////    SET UP COLOR SCALE   /////
        ///////////////////////////////////
        tileColorScale
            .domain(d3.extent(data, colorAccessor));

        ////////////////////////////////////
        /////    ADD CHART ELEMENTS    /////
        ////////////////////////////////////
        annotationGap = tileChartDimensions.boundedWidth * 0.5;
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
        
        // Add horizontal black line to differentiate years
        tileChartBounds.select(".annotations")
            .append("line")
                .attr("x1", 0)
                .attr("x2", chartDimensions.boundedWidth)
                .attr("y1", yScale(370))
                .attr("y2", yScale(370))
                .style("stroke", "#000000")
                .style("stroke-width", 1)
                .style("stroke-dasharray", ("2, 5"));

        // Add year labels
        tileChartBounds.select(".annotations")
            .append("text")
                .attr("class", "axis-title")
                .attr("y", yScale(365))
                .attr("x", annotationGap / 2)
                .attr("text-anchor", "middle")
                .attr("dominant-baseline", "ideographic")
                .text("2016")

        tileChartBounds.select(".annotations")
            .append("text")
                .attr("class", "axis-title")
                .attr("y", yScale(375))
                .attr("x", annotationGap / 2)
                .attr("text-anchor", "middle")
                .attr("dominant-baseline", "hanging")
                .text("2015")
    }

    function addTileLegend() {
        // build list of possible total particle counts
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
            .attr("x", tileChartDimensions.boundedWidth / 2)
            .attr("y", -tileChartDimensions.margin.top)
            .attr("dx", 0)
            .attr("dy", 0)
            .attr("text-anchor", "middle")
            .attr("dominant-baseline", "hanging")
            .attr("text-width", tileChartDimensions.boundedWidth)
            .text("Particulate count")
            .call(d => wrap(d, {shift: false}))

        // append legend rectangle
        const rectWidth = tileChartDimensions.boundedWidth / 2;
        const rectHeight = mobileView ? tileChartDimensions.margin.top / 8 : tileChartDimensions.margin.top / 6;
        const rectX = tileChartDimensions.boundedWidth / 2 - rectWidth / 2;
        legendGroup.append("rect")
            .attr("class", "c1p2 matrixLegend")
            .attr("width", rectWidth)
            .attr("height", rectHeight)
            .attr("fill", "url(#gradient-particles)")
            .attr("x", rectX)
            .attr("y", -tileChartDimensions.margin.top / 2 - rectHeight / 2)

        // append legend text
        const xBuffer = 5;
        legendGroup.append("text")
            .attr("class", "axis-subtitle")
            .attr("text-anchor", "end")
            .attr("dominant-baseline", "central")
            .attr("x", rectX - xBuffer)
            .attr("y", -tileChartDimensions.margin.top / 2)
            .text("low")

        legendGroup.append("text")
            .attr("class", "axis-subtitle")
            .attr("text-anchor", "start")
            .attr("dominant-baseline", "central")
            .attr("x", rectX + rectWidth + xBuffer)
            .attr("y", -tileChartDimensions.margin.top / 2)
            .text("high")
        
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

        barChartBounds.append("text")
            .attr("class", "data-notation")
            .attr("x", 0)
            .attr("y", barYScale('400') + barYScale.bandwidth() / 2)
            .attr("dominant-baseline", "central")
            .text("missing data")
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
        // const interItemSpacing = mobileView ? 15 : 10;
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
            .attr("y", -barChartDimensions.margin.top / 1.75 - legendRectSize / 2.5)
            .attr("width", legendRectSize)
            .attr("height", legendRectSize)
            .style("fill", d => barColorScale(d))
        
        // Add text for each group
        legendGroup.append("text")
            .attr("class", "legend-text")
            .attr("x", legendRectSize + intraItemSpacing) // put text to the right of the rectangle
            .attr("y", -barChartDimensions.margin.top / 1.75)
            .attr("text-anchor", "start") // left-align text
            .attr("dominant-baseline", "central")
            .text(d => d.toLowerCase());

        // Position legend groups
        // https://stackoverflow.com/questions/20224611/d3-position-text-element-dependent-on-length-of-element-before
        // const xBuffer = 6; // set xBuffer for use in mobile row x translations
        legendGroup
            .attr("transform", (d, i) => {
            // Compute total width of preceeding legend items, with spacing
            // Start with width of legend title
            // const d3.select('#legend-title')._groups[0][0].getBBox().width + interItemSpacing;
            
            // Begin right of legend
            const cumulativeWidth = 0;
            // let cumulativeWidth = titleWidth;
            // if (!mobileView) {
            //     // row 1: items 0 and 1
            //     if (i < 2) {
            //         for (let j = 0; j < i; j++) {
            //             cumulativeWidth = cumulativeWidth + legendGroup._groups[0][j].getBBox().width + interItemSpacing;
            //         }
            //     }
            //     // row 2: item 2
            //     else if (i < 3) {
            //         for (let j = 2; j < i; j++) {
            //             cumulativeWidth = cumulativeWidth + legendGroup._groups[0][j].getBBox().width + interItemSpacing;
            //         }
            //     }
            // } else {
                // if (i < 1) {
                //     for (let j = 0; j < i; j++) {
                //         cumulativeWidth = cumulativeWidth + legendGroup._groups[0][j].getBBox().width + interItemSpacing;
                //     }
                // }
                // // row 2: item 2
                // else if (i < 2) {
                //     for (let j = 2; j < i; j++) {
                //         cumulativeWidth = cumulativeWidth + legendGroup._groups[0][j].getBBox().width + interItemSpacing;
                //     }
                // }
                // // row 3: item 3 
                // else if (i === 3) {
                //     for (let j = 2; j < i; j++) {
                //         cumulativeWidth = cumulativeWidth + legendGroup._groups[0][j].getBBox().width + interItemSpacing;
                //     }
                // }
            // }

            // Determine x and y translation
            // set y translation for each row
            let rowHeight = window.innerHeight < 600 ? legendRectSize * 4.5 : legendRectSize * 3;
            rowHeight =  mobileView ? legendRectSize * 4.5 : rowHeight;
            const yTranslation = rowHeight * i;
            // let yTranslation = 0;
            // if (!mobileView) {
            //     if (i < 2) {
            //         yTranslation = 0;
            //     } else if (i < 5) {
            //         yTranslation = window.innerHeight < 770 ? legendRectSize * 4 : legendRectSize * 2;
            //         cumulativeWidth = cumulativeWidth - titleWidth;
            //     } else {
            //         yTranslation = legendRectSize * 5.75
            //         cumulativeWidth = xBuffer; // for last item just translate by xBuffer
            //     } 
            // } else {
                // if (i < 1) {
                //     yTranslation = 0;
                // } else if (i < 2) {
                //     yTranslation = window.innerHeight < 770 ? legendRectSize * 4 : legendRectSize * 2;
                //     cumulativeWidth = cumulativeWidth - titleWidth;
                // } else {
                //     yTranslation = window.innerHeight < 770 ? legendRectSize * 8 : legendRectSize * 4;
                //     cumulativeWidth = 0; // for last item just translate by xBuffer
                // } 
            // }

            // translate each group by that width and height
            return "translate(" + cumulativeWidth + "," + yTranslation + ")"
        })
    }

    function drawScatterChart(data, type) {
        //////////////////////////////
        /////    PROCESS DATA    /////
        //////////////////////////////
        const filteredData = data.filter(d => d.vegetation_type == type)

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
            .domain([... new Set(filteredData.map(d => xAccessor(d)))]);
        drawXAxis(
            {
                axis: scatterXAxis, 
                axisScale: scatterXScale, 
                chartDims: scatterChartDimensions
            }, 
            {
                axisPosition: scatterXAxisPosition, 
                axisTitle: 'Burned vegetation type',
                wrapTitle: true,
                tickType: 'string',
                keepDomain: false,
                keepLabels: false
            }
        )

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
                .data(filteredData, d => d.vegetation_type) // bind data
                .join(
                    enter => enter
                    .append("circle")
                        .attr("class", "point")
                        .attr("id", d => 'point-' + identifierAccessor(d))
                        .attr("cx", d => scatterXScale(xAccessor(d)) + scatterXScale.bandwidth()/2)
                        .attr("cy", d => barYScale(yAccessor(d)) + barYScale.bandwidth()/2)
                        .attr("r", barYScale.bandwidth() / 2 * 0.95)
                        .style("fill", d => scatterColorScale(colorAccessor(d)))
                        .style("fill-opacity", 0)
                        .transition()
                        .duration(transitionLength)
                        .style("fill-opacity", 1),

                    null, // no update function

                    exit => {
                        exit
                        .transition()
                        .duration(transitionLength)
                        .style("fill-opacity", 0)
                        .remove();
                    }
                );

        const baseWidth = mobileView ? scatterChartTranslateX3 - scatterChartDimensions.boundedWidth + (barYScale.bandwidth() / 2 * 0.95 * 2) : scatterChartTranslateX3 - scatterChartDimensions.boundedWidth / 2 + (barYScale.bandwidth() / 2 * 0.95 * 2);
        scatterChartBounds.select('.points') // selects our group we set up to hold chart elements
            .selectAll(".rect") // empty selection
                .data(filteredData, d => d.vegetation_type) // bind data
                .join(
                    enter => enter
                    .append("rect")
                        .attr("class", "rect")
                        .attr("id", d => 'point-' + identifierAccessor(d))
                        .attr("x", d => scatterXScale(xAccessor(d)) + scatterXScale.bandwidth()/2)
                        .attr("y", d => barYScale(yAccessor(d)))
                        .attr("height", barYScale.bandwidth())
                        .attr("width", d => scatterXScale(xAccessor(d)) + baseWidth) //(chartDimensions.boundedWidth - scatterChartDimensions.boundedWidth - chartGap)
                        .attr("transform", d => "translate(" + - (baseWidth + scatterXScale(xAccessor(d))) + ", 0)")
                        .style("fill", d => scatterColorScale(colorAccessor(d)))
                        .style("opacity", 0)
                        .transition()
                        .duration(transitionLength)
                        .style("opacity", 0.5),

                    null, // no update function

                    exit => {
                        exit
                        .transition()
                        .duration(transitionLength)
                        .style("fill-opacity", 0)
                        .remove();
                    }
                );
    }

    function addScatterLegend() {
        const scatterLegendGroup = scatterChartBounds.append("g")
            .attr("id", "bar-chart-legend")

        // Add legend title
         // add axis title
        // scatterLegendGroup.append("text")
        //     .attr("id", "legend-title")
        //     .attr("class", "axis-title")
        //     .attr("x", scatterChartDimensions.boundedWidth / 2)
        //     .attr("y", -scatterChartDimensions.margin.top + 10)
        //     .attr("dx", 0)
        //     .attr("dy", 0)
        //     .attr("text-anchor", "middle")
        //     .attr("dominant-baseline", "hanging")
        //     .attr("text-width", scatterChartDimensions.boundedWidth)
        //     .text('Burned vegetation type')
        //     .call(d => wrap(d, {shift: false}))

        const legendPointSize = barYScale.bandwidth() / 2 * 0.95;
        // const interItemSpacing = mobileView ? 15 : 10;
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

        // Add points for each group
        legendGroups.append("circle")
            .attr("class", "legend-point")
            .attr("cx", 0)
            .attr("cy", -scatterChartDimensions.margin.top / 2 + legendPointSize / 1.5)
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
        legendGroups
            .attr("transform", (d, i) => {
                // Compute total width of preceeding legend items, with spacing
                // Start with width of legend title
                // const titleWidth = d3.select('#legend-title')._groups[0][0].getBBox().width + interItemSpacing;
                
                // Begin right of legend
                const cumulativeWidth = 0;
                // let cumulativeWidth = titleWidth;
                // if (!mobileView) {
                //     // row 1: items 0 and 1
                //     if (i < 2) {
                //         for (let j = 0; j < i; j++) {
                //             cumulativeWidth = cumulativeWidth + legendGroups._groups[0][j].getBBox().width + interItemSpacing;
                //         }
                //     }
                //     // row 2: item 2
                //     else if (i < 3) {
                //         for (let j = 2; j < i; j++) {
                //             cumulativeWidth = cumulativeWidth + legendGroups._groups[0][j].getBBox().width + interItemSpacing;
                //         }
                //     }
                // } else {
                    // if (i < 1) {
                    //     for (let j = 0; j < i; j++) {
                    //         cumulativeWidth = cumulativeWidth + legendGroup._groups[0][j].getBBox().width + interItemSpacing;
                    //     }
                    // }
                    // // row 2: item 2
                    // else if (i < 2) {
                    //     for (let j = 2; j < i; j++) {
                    //         cumulativeWidth = cumulativeWidth + legendGroup._groups[0][j].getBBox().width + interItemSpacing;
                    //     }
                    // }
                    // // row 3: item 3 
                    // else if (i === 3) {
                    //     for (let j = 2; j < i; j++) {
                    //         cumulativeWidth = cumulativeWidth + legendGroup._groups[0][j].getBBox().width + interItemSpacing;
                    //     }
                    // }
                // }

                // Determine x and y translation
                // set y translation for each row               
                let rowHeight = window.innerHeight < 600 ? barYScale.bandwidth() * 4.5 : barYScale.bandwidth() * 3;                const yTranslation = rowHeight * i;
                rowHeight =  mobileView ? barYScale.bandwidth() * 4.5 : rowHeight;
                // let yTranslation = 0;
                // if (!mobileView) {
                //     if (i < 2) {
                //         yTranslation = 0;
                //     } else if (i < 3) {
                //         yTranslation = window.innerHeight < 770 ? barYScale.bandwidth() * 4 : barYScale.bandwidth() * 2;
                //         cumulativeWidth = cumulativeWidth - titleWidth;
                //     }
                // } else {
                    // if (i < 1) {
                    //     yTranslation = 0;
                    // } else if (i < 2) {
                    //     yTranslation = barYScale.bandwidth() * 4;
                    //     cumulativeWidth = cumulativeWidth - titleWidth;
                    // } else {
                    //     yTranslation = barYScale.bandwidth() * 8
                    //     cumulativeWidth = 0; // for last item just translate by xBuffer
                    // }
                // }

                // translate each group by that width and height
                return "translate(" + cumulativeWidth + "," + yTranslation + ")"
            })
    }

    function showChart(hiddenChart, translateX) {
        hiddenChart
            .attr("opacity", 0)
            .attr("visibility", "visible")
            // .attr("display", "auto")
            .transition()
            .duration(transitionLength)
            .attr("opacity", 1)
            .style("transform", `translate(${
                translateX
            }px, 0px)`);
    }

    function hideChart(visibleChart) {
        visibleChart
            .transition()
            .duration(transitionLength)
            .attr("opacity", 0)
        setTimeout(function() {
            visibleChart
                .attr("visibility", "hidden");
        }, transitionLength)
    }

    function moveChart(visibleChart, translateX) {
        visibleChart
            .transition()
            .duration(transitionLength)
            .style("transform", `translate(${
                translateX
            }px, 0px)`);
    }

    function updateChartView(index) {
        const barChartHidden = barChartWrapper.node().attributes.visibility.value == 'hidden';
        const scatterChartHidden = scatterChartWrapper.node().attributes.visibility.value == 'hidden';
        if (index == 1) {
            maskingRect
                .transition()
                .duration(transitionLength)
                .delay(transitionLength / 4)
                .style("transform", `translate(${
                    tileChartTranslateX1 + tileChartDimensions.width + chartGap
                }px, 0px)`);
            moveChart(tileChartWrapper, tileChartTranslateX1)
            if (!barChartHidden) hideChart(barChartWrapper)            
        } else if (index == 2) {
            maskingRect
                .style("transform", `translate(${
                    tileChartTranslateX2 + tileChartDimensions.width + barChartDimensions.width + chartGap * 2.25
                }px, 0px)`)
                .transition()
                .duration(transitionLength)
                .delay(transitionLength / 4)
                .style("opacity", "1");
            moveChart(tileChartWrapper, tileChartTranslateX2)
            if (barChartHidden)  {
                showChart(barChartWrapper, barChartTranslateX2)
            } else {
                moveChart(barChartWrapper, barChartTranslateX2)
            }
            if (!scatterChartHidden) hideChart(scatterChartWrapper) 
        } else if (index == 3) {
            drawScatterChart(scatterData.value, 'softwood')
            maskingRect
                .style("opacity", "0")
            moveChart(tileChartWrapper, tileChartTranslateX3)
            moveChart(barChartWrapper, barChartTranslateX3)
            if (scatterChartHidden)  {
                showChart(scatterChartWrapper, scatterChartTranslateX3)
            }
        } else if (index == 4) {
            maskingRect
                .style("opacity", "0")
            drawScatterChart(scatterData.value, 'hardwood')
        }
    }

    // https://gist.github.com/mbostock/7555321
    function wrap(text, {
        shift = false
    }) {
        text.each(function() {
            var text = d3.select(this),
            words = text.text().split(/\s|-+/).reverse(),
            word,
            line = [],
            lineNumber = 0,
            lineHeight = 1.1, // ems
            width = text.attr("text-width"),
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
            // Likely only want to shift if dominant-baseline = central
            if (lineNumber > 0 && shift) {
                const startDy = -(lineNumber * (lineHeight / 2));
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
        max-width: 900px;
        grid-template-columns: 10% calc(80% - 4rem) 10%;
        grid-template-rows: auto max-content;
        grid-template-areas:
            "prev-upper text next-upper"
            "prev-lower chart next-lower";
        margin: 2rem auto 4rem auto;
        column-gap: 2rem;
        row-gap: 3rem;
        @media only screen and (max-width: 600px) {
            width: 90vw;
            grid-template-rows: auto max-content;
            grid-template-areas:
                "chart chart chart"
                "prev-upper text next-upper";
        }
    }
    #chart-container {
        grid-area: chart;
        width: 100%;
    }
    #aerosol-text-container {
        grid-area: text;
        height: 19vh;
        align-content: center;
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
            height: 3.5rem;
            width: 3.5rem;
            align-self: start;
        }
    }
    #aerosol-prev-upper {
        grid-area: prev-upper;
        justify-self: end;
    }
    #aerosol-next-upper {
        grid-area: next-upper;
        justify-self: start;
    }
    #aerosol-prev-lower {
        grid-area: prev-lower;
        justify-self: end;
    }
    #aerosol-next-lower {
        grid-area: next-lower;
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
    #softwoods-tooltip {
        margin-left: -145px;
    }
    #hardwoods-tooltip {
        margin-left: -145px;
    }
    /* css for elements added/classed w/ d3 */
    #masking-rect {
        fill: var(--color-background);
    }
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
            font-size: 1.4rem;
        }
    }
    .axis-subtitle {
        font-size: 1.8rem;
        font-family: var(--default-font);
        font-weight: 300;
        fill: var(--color-text);
        user-select: none;
        @media only screen and (max-width: 600px) {
            font-size: 1.4rem;
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
    .tick-line {
        stroke: var(--dark-grey);
        stroke-width: 0.5;
        stroke-dasharray: 1px, 2px;
    }
    .data-notation {
        font-size: 1rem;
        font-family: var(--default-font);
        user-select: none;
        @media only screen and (max-width: 600px) {
            font-size: 0.8rem;
        }
        @media screen and (max-height: 770px) { 
            font-size: 0.8rem;
        }
    }
</style>