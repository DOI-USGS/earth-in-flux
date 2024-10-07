<template>
    <!---VizSection-->
    <VizSection
        :figures="true"
        :fig-caption="false"
    >
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
    const dataFile = 'fish_as_food_country_climate.csv'
    const data = ref();
    const chart = ref(null);
    const chartTitle = 'Title of chart';
    let chartDimensions;
    let chartBounds;
    let xScale;
    let xAxis;
    let yScale;
    let yAxis;
    let rScale;
    const hex1 = "#6F4E37" // lower end
    const hex2 =  "#5CB270" // upper end
    const colors = {"1. High income: OECD": mixHexColors(hex1,hex2,1.0),
                    "2. High income: nonOECD": mixHexColors(hex1,hex2,0.75),
                    "3. Upper middle income": mixHexColors(hex1,hex2,0.50),
                    "4. Lower middle income": mixHexColors(hex1,hex2,0.25),
                    "5. Low income": mixHexColors(hex1,hex2,0.0)}
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
                    height: window.innerHeight * 0.6,
                    margin: 10,
                    marginBottom: 50,
                    marginLeft: 70});

                // draw chart
                drawChart(data.value, 'Europe');
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
                d.population = +d.population;
                d.n_fishers = +d.n_fishers;
                d.total_rec_harvest_kg = +d.total_rec_harvest_kg;
                d.total_consumable_harv_kg = +d.total_consumable_harv_kg;
                d.MCDM_VUL_2075_45 = +d.MCDM_VUL_2075_45;
                d.consum_kg_person = +d.consum_kg_person;
                d.consum_kg_fisher = +d.consum_kg_fisher;
                return d;
            });
            return data;
        } catch (error) {
            console.error(`Error loading data from ${fileName}`, error);
            return [];
        }
    }

    function mixHexColors(hex1, hex2, ratio = 0.5) {
        // Helper function to convert hex to RGB
        function hexToRgb(hex) {
            hex = hex.replace(/^#/, '');
            let bigint = parseInt(hex, 16);
            let r = (bigint >> 16) & 255;
            let g = (bigint >> 8) & 255;
            let b = bigint & 255;
            return { r, g, b };
        }

        // Helper function to convert RGB to hex
        function rgbToHex(r, g, b) {
            return `#${((1 << 24) + (r << 16) + (g << 8) + b).toString(16).slice(1).toUpperCase()}`;
        }

        // Convert hex colors to RGB
        const rgb1 = hexToRgb(hex1);
        const rgb2 = hexToRgb(hex2);

        // Mix the colors based on the specified ratio
        const mixedR = Math.round(rgb1.r * (1 - ratio) + rgb2.r * ratio);
        const mixedG = Math.round(rgb1.g * (1 - ratio) + rgb2.g * ratio);
        const mixedB = Math.round(rgb1.b * (1 - ratio) + rgb2.b * ratio);

        // Convert mixed RGB back to hex
        return rgbToHex(mixedR, mixedG, mixedB);
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
        initXScale()
        initYScale()

        // Initialize axes
        initXAxis()
        initYAxis()

        // Add groups for visual elements
        chartBounds.append("g")
            .attr("class", "circles");
    }

    function initXScale() {
        // scale for the x axis (domain set in `drawChart()`)
        xScale = d3.scaleLinear()
            .range([0, chartDimensions.boundedWidth])
    }

    function initXAxis() {
        // add group for x axis
        xAxis = chartBounds.append("g")
            .attr("id", "x-axis")
            .attr("class", "axis")
            .attr("transform", `translate(0,${chartDimensions.boundedHeight})`)
            .attr("aria-hidden", true); // hide from screen reader
    }

    function initYScale() {
        // scale for y axis (domain set in `drawChart()`)
        yScale = d3.scaleLog()
            .range([chartDimensions.boundedHeight, 0]);
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
        textAngle = 0,
        keepDomain = true,
    }) {
        // generate axis
        // if numeric ticks, include specification of format
        if (tickType == "numeric") {
            axis
                .call(d3[axisFxn](axisScale).tickSize(tickSize).tickPadding(tickPadding).tickFormat(d3.format(tickFormat)));
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

    function drawXAxis({
        axisTitle = '',
        titleX = chartDimensions.boundedWidth / 2,
        titleY = chartDimensions.margin.bottom,
        titleTextAnchor = "middle",
        titleBaseline = "text-after-edge",
        titleAngle = 0,
        tickSize = 0,
        tickPadding = 5,
        tickType = 'numeric',
        tickFormat = ".2f",
        textAngle = 0, 
        keepDomain = true,
    }) {
        drawAxis({
            axis: xAxis,
            axisScale: xScale,
            axisFxn: 'axisBottom'
        }, {
            axisTitle: axisTitle,
            titleX: titleX,
            titleY: titleY,
            titleTextAnchor: titleTextAnchor,
            titleBaseline: titleBaseline,
            titleAngle: titleAngle,tickSize: tickSize,
            tickPadding: tickPadding,
            tickType: tickType,
            tickFormat: tickFormat,
            textAngle: textAngle,
            keepDomain: keepDomain,
        })
    }

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
            keepDomain: keepDomain,
        })
    }

    function initRScale() {
        rScale = d3.scaleSqrt()
            .range([chartDimensions.boundedWidth*0.005,chartDimensions.boundedWidth*0.05]); // from 0.5% to 5% of chart width
    }

    function initColorScale(data) {
        colorScale = d3.scaleOrdinal()
            .domain(data)
            .range(data.map(item => colors[item]));
    }

    function drawChart(data, continent) {
        //////////////////////////////
        /////    PROCESS DATA    /////
        //////////////////////////////
        const chartData = data //.filter(d => d.continent === continent); // May filter later

        ///////////////////////////////////////////
        /////    SET UP ACCESSOR FUNCTIONS    /////
        ///////////////////////////////////////////
        const xAccessor = d => d.MCDM_VUL_2075_45;
        const yAccessor = d => d.consum_kg_person;
        const rAccessor = d => d.population;
        const colorAccessor = d => d.income_grp;
        const identifierAccessor = d => d.admin.replace(/ /g,"_");

        ///////////////////////////////////////////
        /////    FINISH SETTING UP X SCALE    /////
        ///////////////////////////////////////////
        // set domain for xScale, based on data
        xScale
            .domain(d3.extent(chartData, xAccessor));
        drawXAxis({axisTitle: 'Climate vulnerability'})
        
        ///////////////////////////////////////////
        /////    FINISH SETTING UP Y SCALE    /////
        ///////////////////////////////////////////
        // set domain for yScale
        yScale
            .domain(d3.extent(chartData, yAccessor));
        drawYAxis({axisTitle: 'Per capita consumption, in kilograms', tickFormat:".0e"})

        ///////////////////////////////////////////
        /////    FINISH SETTING UP R SCALE    /////
        ///////////////////////////////////////////
        // set domain for rScale
        initRScale()
        rScale
            .domain(d3.extent(chartData, rAccessor));
            
        // ///////////////////////////////////
        // /////    SET UP COLOR SCALE   /////
        // ///////////////////////////////////
        const colorCategories = [... new Set(data.map(colorAccessor))];
        initColorScale(colorCategories)

        ////////////////////////////////////
        /////    ADD CHART ELEMENTS    /////
        ////////////////////////////////////
        // draw chart
        chartBounds.select('.circles') // selects our group we set up to hold chart elements
            .selectAll(".circle") // empty selection
                .data(chartData) // bind data
                .enter() // instantiate chart element for each element of data
                .append("circle") // append a circle for each element
                    .attr("class", d => "circle " + d.continent)
                    .attr("id", d => 'circle-' + identifierAccessor(d))
                    .attr("cx", d => xScale(xAccessor(d)))
                    .attr("cy", d => yScale(yAccessor(d)))
                    .attr("r", d => rScale(rAccessor(d)))
                    .attr('fill', d => colorScale(colorAccessor(d)));

    }
</script>

<style scoped lang="scss">
</style>

<style lang="scss">
/* css for elements added and classed w/ d3 */
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
</style>