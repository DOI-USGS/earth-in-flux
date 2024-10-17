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
        <template #aboveExplanation>
            <p v-html="text.paragraph1" />
            <p v-html="text.paragraph2" />
            <p v-html="text.paragraph3" />
            <p v-html="text.paragraph4" />
        </template>
        <template #figures>
            <div id="chart-container" class="maxWidth" ref="chart"></div>
        </template>
        <template #figureCaption>
        </template>
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
    const dataFile = 'fii_core4data.csv'
    const data = ref();
    const chart = ref(null);
    const chartTitle = 'Title of chart';
    let chartDimensions;
    let chartBounds;
    // let xScale;
    // let xAxis;
    let yScale;
    let yAxis;
    let colorScale;

    // Behavior on mounted (functions called here)
    // Load data and then make chart
    onMounted(async () => {
        try {
            await loadDatasets();

            if (data.value.length > 0) {
                // initialize chart elements
                initChart({
                    width: chart.value.offsetWidth,
                    height: window.innerHeight * 0.8,
                    margin: 10,
                    marginBottom: 50,
                    marginLeft: 100});

                // draw chart
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
            const data = await d3.csv(publicPath + fileName, d => {
                d.depth_cm = +d.depth_cm;
                d.total = +d.total;
                return d;
            });
            return data;
        } catch (error) {
            console.error(`Error loading data from ${fileName}`, error);
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

        // Initialize scales
        // initXScale()
        initYScale()
        initColorScale()

        // Initialize axes
        // initXAxis()
        initYAxis()

        // Add groups for visual elements
        chartBounds.append("g")
            .attr("class", "rects");
        chartBounds.append("g")
            .attr("class", "annotations");
    }

    // function initXScale() {
    //     // scale for x axis (domain set in `drawChart()`)
    //     xScale = d3.scaleLinear()
    //         .range([0, chartDimensions.boundedWidth]);
    // }

    // function initXAxis() {
    //     // add group for x axis
    //     xAxis = chartBounds.append("g")
    //         .attr("id", "x-axis")
    //         .attr("class", "axis")
    //         .attr("transform", `translate(0,${chartDimensions.boundedHeight})`)
    //         .attr("aria-hidden", true); // hide from screen reader
    // }

    function initYScale() {
        // scale for the y axis (domain set in `drawChart()`)
        yScale = d3.scaleLinear()
            .range([0, chartDimensions.boundedHeight]);
    }

    function initYAxis() {
        // add group for y axis
        yAxis = chartBounds.append("g")
            .attr("id", "y-axis")
            .attr("class", "axis")
            .attr("aria-hidden", true);
    }

    function drawAxis({
        axis,
        axisScale,
        axisFxn
    }, {
        axisTitle = '',
        titleX = -chartDimensions.boundedHeight / 2,
        titleY = -chartDimensions.margin.left,
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
    //     titleX = chartDimensions.boundedWidth / 2,
    //     titleY = chartDimensions.margin.bottom,
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
        axisTitle = '',
        titleX = -chartDimensions.boundedHeight / 2,
        titleY = -chartDimensions.margin.left,
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
            axis: yAxis,
            axisScale: yScale,
            axisFxn: 'axisLeft'
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

    function initColorScale() {
        colorScale = d3.scaleSequential()           
            .interpolator(d3.interpolateGreys);
    }

    function drawChart(data) {
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
        drawYAxis({tickFormat: ".0f", customSuffix: 'cm', tickSize: 3, keepDomain: false})
        
        ///////////////////////////////////////////
        /////    FINISH SETTING UP X SCALE    /////
        ///////////////////////////////////////////
        // // set domain for xScale
        // xScale
        //     .domain([0, d3.max(data, xAccessor)]);
        // drawXAxis({axisTitle: 'Climate vulnerability'})

        ///////////////////////////////////
        /////    SET UP COLOR SCALE   /////
        ///////////////////////////////////
        colorScale
            .domain(d3.extent(data, colorAccessor));

        ////////////////////////////////////
        /////    ADD CHART ELEMENTS    /////
        ////////////////////////////////////
        // draw chart
        chartBounds.select('.rects') // selects our group we set up to hold chart elements
            .selectAll(".rect") // empty selection
                .data(data) // bind data
                .enter() // instantiate chart element for each element of data
                .append("rect") // append a rectangle for each element
                    .attr("class", "rect")
                    .attr("id", d => 'bar-' + identifierAccessor(d))
                    .attr("x", 40)
                    .attr("y", d => yScale(yAccessor(d)))
                    .attr("height", chartDimensions.boundedHeight / data.length)
                    .attr("width", 50)
                    .style("fill", d => colorScale(colorAccessor(d)));
        // draw year bands
        chartBounds.select(".annotations")
            .append("rect")
                .attr("class", "year-bands")
                .attr("x", 25)
                .attr("y", 0)
                .attr("height", yScale(372))
                .attr("width", 3)

        chartBounds.select(".annotations")
            .append("rect")
                .attr("class", "year-bands")
                .attr("x", 25)
                .attr("y", yScale(378))
                .attr("height", chartDimensions.boundedHeight - yScale(378))
                .attr("width", 3)

        chartBounds.select(".annotations")
            .append("text")
                .attr("x", - yScale(372) / 2)
                .attr("y", 20)
                .attr("transform", "rotate(-90)")
                .attr("text-anchor", "middle")
                .text("2016 accumulation")
        
        chartBounds.select(".annotations")
            .append("text")
                .attr("x", - yScale(378) - ((chartDimensions.boundedHeight - yScale(378)) / 2))
                .attr("y", 20)
                .attr("transform", "rotate(-90)")
                .attr("text-anchor", "middle")
                .text("2015 accumulation")
    }
</script>

<style scoped lang="scss">
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