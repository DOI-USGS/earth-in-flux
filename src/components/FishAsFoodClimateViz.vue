<template>
  <!---VizSection-->
  <VizSection
    :figures="true"
    :fig-caption="false"
  >
    <!-- EXPLANATION -->
    <template #aboveExplanation>
      <p v-html="text.paragraph1" />
    </template>
    <template #figures>
      <div
        id="grid-container"
      >
        <CountryInfoBox :active-country="activeCountry" />
        <div
          class="axis-title-container"
          aria-hidden="true"
        >
          <p class="axis-title">
            {{ yAxisLabel }}
          </p>
        </div>
        <div
          id="chart-container"
          ref="chart"
          class="maxWidth"
        />
      </div>
    </template>
    <!-- FIGURE CAPTION -->
    <template #figureCaption />
  </VizSection>
</template>

<script setup>
    import { computed, onMounted, ref } from "vue";    
    import { isMobile } from 'mobile-device-detect';
    import * as d3 from 'd3';
    import VizSection from '@/components/VizSection.vue';
    import CountryInfoBox from '@/components/CountryInfoBox.vue';

    // define props
    const props = defineProps({
      text: { 
        type: Object,
        default: () => ({}),
      }
    })

    // global variables
    const publicPath = import.meta.env.BASE_URL;
    const mobileView = isMobile;
    const defaultInfo = props.text.defaultInfo;
    const activeCountry = ref(props.text.defaultInfo);
    const data = ref();
    const datasetsConfig = [
        { file: 'fish_as_food_country_climate.csv', ref: data, type: 'csv', numericFields: ['population','participation_rate','n_fishers', 'total_rec_harv_kg','total_consumable_harv_kg', 'MCDM_VUL_2075_45', 'consum_kg_fisher', 'warm', 'cool', 'cold', 'unknown']}
    ]
    const chart = ref(null);
    const currentVal = ref(null);
    const currentCountry = ref(null);
    const countrySelected = ref(false);
    const chartTitle = 'Scatterplot of inland fish consumption by recreational fishers vs. average climate vulnerability of harvested species. There is a circle for each country, and the size of each circle is proportional to the number of recreational fishers in each country. Hit the spacebar or enter while on this element to return to the default view after interacting with the chart.';
    let chartDimensions;
    let chartSVG;
    let chartBounds;
    let xScale;
    let xAxis;
    let yScale;
    let yAxis;
    let rScale;
    const rPropMin = 0.01
    const rPropMax = 0.05
    const colors = {
        Africa: "var(--mid-teal)",
        Oceania: "var(--mid-purple)",
        "North America": "var(--mid-red)",
        Europe: "var(--mid-brown)",
        Asia: "var(--dark-blue)",
        "South America": "var(--grey-blue)"
    }

    const yAxisLabel = computed(() => {
        const defaultLabel = 'Each year, each inland recreational fisher consumes'
        return currentCountry.value ? `Each year, each inland recreational fisher in ${currentCountry.value} consumes ${currentVal.value} of fish` : `${defaultLabel}...`
    })

    // Behavior on mounted (functions called here)
    // Load data and then make chart
    onMounted(async () => {
        try {
            await loadDatasets(datasetsConfig);

            if (data.value.length > 0) {
                // initialize chart elements
                initChart({
                    width: chart.value.offsetWidth,
                    height: mobileView ? window.innerHeight * 0.8 : window.innerHeight * 0.7,
                    margin: 10,
                    marginBottom: mobileView ? 90: 60,
                    marginLeft:5,
                    marginTop: mobileView ? 30 : 25
                });

                // draw chart
                drawChart(data.value);
            } else {
                console.error('Error loading data');
            }
        } catch (error) {
            console.error('Error during component mounting', error);
        }
    });

    async function loadDatasets(configs) {
        for (const { file, ref, type, numericFields} of configs) {
        try {
            ref.value = await loadData(file, type, numericFields);
        } catch (error) {
            console.error(`Error loading ${file}`, error);
        }
        }
    }

    async function loadData(dataFile, dataType, dataNumericFields) {
        try {
          let data;
          if (dataType == 'csv') {
            data = await d3.csv(publicPath + dataFile, d => {
            if (dataNumericFields) {
                dataNumericFields.forEach(numericField => {
                d[numericField] = +d[numericField]
                });
            }
            return d;
            });
          } else if (dataType == 'json') {
            data = await d3.json(publicPath + dataFile);
          } else {
            console.error(`Data type ${dataType} is not supported. Data type must be 'csv' or 'json'`)
          }

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
        chartSVG = d3.select("#chart-container")
            .append("svg")
                .attr("viewBox", [0, 0, (chartDimensions.width), (chartDimensions.height)].join(' '))
                .attr("width", "100%")
                .attr("height", "100%")
                .attr("id", "wrapper");

        // assign aria-label for accessibility
        chartSVG
            .attr("aria-label", chartTitle);

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
            .attr("class", "circles")
            .attr("role", "list")
            .attr("aria-label", "circles for individual countries");
    }

    function initXScale() {
        // scale for the x axis (domain set in `drawChart()`)
        xScale = d3.scaleLinear()
            .range([0, chartDimensions.boundedWidth])
    }

    function initXAxis() {
        // add group for x axis
        xAxis = chartBounds.append("g")
            .attr("aria-hidden", "true")
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
            .attr("aria-hidden", "true")
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
        customFormat = false,
        tickValues = null,
        textAngle = 0,
        keepDomain = true,
        wrapTitle = false,
        wrapWidth = null
    }) {
        // generate axis
        // if numeric ticks, include specification of format
        if (tickType == "numeric" && tickValues && customFormat) {
            axis
                .call(d3[axisFxn](axisScale).tickValues(tickValues).tickSize(tickSize).tickPadding(tickPadding).tickFormat(tickFormat));
        }
        else if (tickType == "numeric" && tickValues && !customFormat) {
            axis
                .call(d3[axisFxn](axisScale).tickValues(tickValues).tickSize(tickSize).tickPadding(tickPadding).tickFormat(d3.format(tickFormat)));
        }
        else if (tickType == "numeric" && customFormat) {
            axis
                .call(d3[axisFxn](axisScale).tickSize(tickSize).tickPadding(tickPadding).tickFormat(tickFormat));
        }
        else if (tickType == "numeric") {
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
            .attr("dx", "0em")
            .attr("dy", wrapTitle ? "-1em" :  "0em")
            .attr("transform", `rotate(${titleAngle})`)
            .attr("text-anchor", titleTextAnchor)
            .attr("dominant-baseline", titleBaseline)
            .attr("role", "presentation")
            .attr("aria-hidden", true)
            .text(axisTitle)
            .attr("text-width", wrapTitle ? wrapWidth : null)
            .call(d => wrapTitle ? wrap(d, {shift: false}) : d);
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
        customFormat = false,
        tickValues = null,
        textAngle = 0, 
        keepDomain = true,
        wrapTitle = false,
        chartWidth = null
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
            customFormat: customFormat,
            tickValues: tickValues,
            textAngle: textAngle,
            keepDomain: keepDomain,
            wrapTitle: wrapTitle,
            wrapWidth: chartWidth
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
        customFormat = false,
        tickValues = null,
        textAngle = 0,
        keepDomain = true,
        wrapTitle = false,
        chartHeight = null
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
            customFormat: customFormat,
            tickValues: tickValues,
            keepDomain: keepDomain,
            wrapTitle: wrapTitle,
            wrapWidth: chartHeight
        })
    }

    function initRScale() {
        rScale = d3.scaleSqrt()
            .range([Math.min(chartDimensions.boundedWidth,chartDimensions.boundedHeight)*rPropMin,Math.min(chartDimensions.boundedWidth,chartDimensions.boundedHeight)*rPropMax]);
    }

    function drawChart(data) {
        //////////////////////////////
        /////    PROCESS DATA    /////
        //////////////////////////////
        const chartData = data.sort((a, b) => b.n_fishers - a.n_fishers) //.filter(d => d.continent === continent); // May filter later

        ///////////////////////////////////////////
        /////    SET UP ACCESSOR FUNCTIONS    /////
        ///////////////////////////////////////////
        const xAccessor = d => d.MCDM_VUL_2075_45;
        const yAccessor = d => d.consum_kg_fisher;
        const rAccessor = d => d.n_fishers;
        const colorAccessor = d => d.continent;
        const identifierAccessor = d => d.admin.replace(/ /g,"_");

        ///////////////////////////////////////////
        /////    FINISH SETTING UP X SCALE    /////
        ///////////////////////////////////////////
        // set domain for xScale, based on data

        // set up x-scale vars
        const xMin = d3.min(chartData, xAccessor);
        const xMax = d3.max(chartData, xAccessor);
        const xRange = xMax - xMin;

        // buffer logic
        // ensures that the buffer is the max radius away from the smallest data point
        const xBuffer = (xRange * rPropMax) / (1 - 2 * rPropMax);

        // domain vars
        let xDomain_min = xMin - xBuffer;
        xDomain_min = mobileView ? xDomain_min * 0.95 : xDomain_min;
        const xDomain_max = xMax + xBuffer;     
        
        // apply domain to scale
        xScale
            .domain([xDomain_min, xDomain_max]);
        drawXAxis({axisTitle: 'Average climate vulnerability of harvested species', tickValues: [], titleBaseline: 'text-after-edge', wrapTitle: true, chartWidth: chartDimensions.boundedWidth})

        // Add arrow
        const arrow_dim = 12;
        chartSVG.append("defs").attr("aria-hidden", "true").append("marker")
            .attr("id", "arrowhead-right") // Assign a unique ID
            .attr("viewBox", "0 -5 10 10") // Define the coordinate system for the marker
            .attr("refX", 10) // Reference point for the end of the arrow
            .attr("refY", 0)
            .attr("markerWidth", arrow_dim)
            .attr("markerHeight", arrow_dim)
            .append("path")
            .attr("d", "M 0,-5 L 10,0 L 0,5"); // Define the arrowhead shape
        chartSVG.append("defs").attr("aria-hidden", "true").append("marker")
            .attr("id", "arrowhead-left") // Assign a unique ID
            .attr("viewBox", "0 -5 10 10") // Define the coordinate system for the marker
            .attr("refX", 0) // Reference point for the end of the arrow
            .attr("refY", 0)
            .attr("markerWidth", arrow_dim)
            .attr("markerHeight", arrow_dim)
            .append("path")
            .attr("d", "M 0,0 L 10,5 L 10,-5"); 
        
        d3.select("#x-axis").select("path.domain")
          .attr("marker-end", "url(#arrowhead-right)") // Reference the defined marker by ID
          .attr("marker-start", "url(#arrowhead-left)"); // Reference the defined marker by ID

        d3.select("#x-axis").append("text")
            .attr("class", "axis-subtitle")
            .attr("x", 0)
            .attr("y", mobileView ? chartDimensions.margin.bottom / 3 : chartDimensions.margin.bottom / 2)
            .attr("dominant-baseline", "text-after-edge")
            .style("text-anchor", "start")
            .text("lower​")

        d3.select("#x-axis").append("text")
            .attr("class", "axis-subtitle")
            .attr("x", chartDimensions.boundedWidth)
            .attr("y", mobileView ? chartDimensions.margin.bottom / 3 : chartDimensions.margin.bottom / 2)
            .attr("dominant-baseline", "text-after-edge")
            .style("text-anchor", "end")
            .text("higher")

        ///////////////////////////////////////////
        /////    FINISH SETTING UP Y SCALE    /////
        ///////////////////////////////////////////
        // set domain for yScale
        const yInnerDomainRange = Math.log10(d3.max(chartData, yAccessor)) - Math.log10(d3.min(chartData, yAccessor));
        const yDomain_min = Math.pow(10,Math.log10(d3.min(chartData, yAccessor)) - yInnerDomainRange * rPropMax / (1.0 - 2.0 * rPropMax)); //ensures that the buffer is the max radius away from the smallest data point
        const yDomain_max = Math.pow(10,Math.log10(d3.max(chartData, yAccessor)) + yInnerDomainRange * rPropMax / (1.0 - 2.0 * rPropMax)); //ensures that the buffer is the max radius away from the larger data point  
        yScale
            .domain([yDomain_min, yDomain_max]);
        const tickFormatter = function(d) {
            let suffix;
            if (d < 1) {
              suffix = "g of fish"              
              return d3.format(".1s")(d).replace("m","") + " " + suffix;
            } else {
              suffix = "kg of fish"              
              return d3.format(".1s")(d) + " " + suffix;
            }            
        }
        drawYAxis({tickFormat: tickFormatter, customFormat: true, tickValues: [0.01, 0.1, 1, 10, 20], tickSize: - chartDimensions.boundedWidth, keepDomain: false})
        // Further customize axis label placement
        d3.select("#y-axis").selectAll("g").selectAll("text")
            .attr("dx", "0.2em")
            .attr("dy", mobileView ? "-0.05em" : "-0.2em")
            .attr("text-anchor", "start")
            .attr("dominant-baseline", "text-after-edge")

        ///////////////////////////////////////////
        /////    FINISH SETTING UP R SCALE    /////
        ///////////////////////////////////////////
        // set domain for rScale
        initRScale()
        rScale
            .domain(d3.extent(chartData, rAccessor));

        ////////////////////////////////////
        /////    ADD CHART ELEMENTS    /////
        ////////////////////////////////////
        // draw chart

        chartBounds.select('.circles') // selects our group we set up to hold chart elements
            .selectAll(".circle") // empty selection
                .data(chartData) // bind data
                .enter() // instantiate chart element for each element of data
                .append("circle") // append a circle for each element
                    .attr("aria-hidden", "true")
                    .attr("class", d => "circle " + d.continent)
                    .attr("id", d => 'circle-' + identifierAccessor(d))
                    .attr("cx", d => xScale(xAccessor(d)))
                    .attr("cy", d => yScale(yAccessor(d)))
                    .attr("r", d => rScale(rAccessor(d)))
                    .attr('fill', d => colors[colorAccessor(d)])
                    .attr('stroke', '#FFFFFF')
                    .attr('stroke-width', '1px')

        chartBounds.select('.circles') // selects our group we set up to hold chart elements
            .selectAll(".circle-mask") // empty selection
                .data(chartData) // bind data
                .enter() // instantiate chart element for each element of data
                .append("circle") // append a circle for each element
                    .attr("tabindex", 0)
                    .attr("role", "listitem")
                    .attr("aria-label", d => {
                        const value = yAccessor(d)
                        const roundedValue = value < 1 ? (value * 1000 < 10 ? `${Math.round(value * 10000) / 10} g` : `${Math.round(value * 1000 / 5) * 5} g`) : `${Math.round(value * 10) / 10} kg`
                        return `Each year, each inland recreational fisher in ${d.admin} consumes ${roundedValue} of fish`
                    })
                    .attr("class", d => "circle-mask " + d.continent)
                    .attr("id", d => 'circle-mask-' + identifierAccessor(d))
                    .attr("cx", d => xScale(xAccessor(d)))
                    .attr("cy", d => yScale(yAccessor(d)))
                    .attr("r", d => rScale(rAccessor(d)))
                    .attr('fill', '#FFFFFF')
                    .attr('stroke', '#FFFFFF')
                    .attr('stroke-width', '1px')
                    .attr('fill-opacity', 0)
                    .attr('stroke-opacity', 0)
                    .on("click", (event, d) => {
                        selectCountry(event, d)
                    })
                    .on("keydown", function(event, d) {
                        if (event.code == 'Enter' | event.code == 'Space') {
                          selectCountry(event, d)
                        }
                    })
                    .on("mouseover", (event, d) => {
                        // If no country is selected, allow mouseover updates to y-axis label
                        if (!countrySelected.value) {
                            populateYAxisLabel(d.admin, yAccessor(d))
                        }
                    })
                    .on("mouseout", () => {
                        // If no country is selected, allow mouseout reset of y-axis label
                        if (!countrySelected.value) {
                            currentCountry.value = null;
                        }
                    })
        
        function selectCountry(event, d) {
            // Note that a country is selected
            countrySelected.value = true;
            // Update the y-axis label
            populateYAxisLabel(d.admin, yAccessor(d))
            
            chartBounds.selectAll(".circle")
                .attr('stroke', '#FFFFFF')
                .attr('stroke-width', '1px')
            chartBounds.selectAll(".circle-mask")
                .attr("fill-opacity", 0.5)
                .attr("stroke-opacity", 0.5)
            // arrange circles in DOM
            for (var i = 0; i < chartData.length; i++){
                let currentCircle = chartBounds.select("#circle-" + identifierAccessor(chartData[i]))
                let currentCircleMask = chartBounds.select("#circle-mask-" + identifierAccessor(chartData[i]))
                currentCircle.raise();
                currentCircleMask.raise();
            };
            chartBounds.select("#circle-" + identifierAccessor(d))
                .attr('stroke', '#000000')
                .attr('stroke-width', '3px')
                .raise()
            chartBounds.select("#circle-mask-" + identifierAccessor(d))
                .attr("fill-opacity", 0)
                .attr("stroke-opacity", 0)
                .raise()
            activeCountry.value = {
                name: d.admin,
                continent: colorAccessor(d).replace(/ /g,"_"),
                participation_rate: Math.round(d.participation_rate),
                n_fishers: d3.format(',')(d.n_fishers),
                consum_harv_kg: d3.format(',')(Math.round(d.total_consumable_harv_kg)),
                mean_vul: Math.round(d.MCDM_VUL_2075_45*100)/100,
                harv_breakdown: [
                    {
                        guild: "warm",
                        percent: Math.round(d.warm*10)/10
                    },
                    {
                        guild: "cool",
                        percent: Math.round(d.cool*10)/10
                    },
                    {
                        guild: "cold",
                        percent: Math.round(d.cold*10)/10
                    },
                    {
                        guild: "unknown",
                        percent: Math.round(d.unknown*10)/10
                    }
                ]
            }
            // Don't propagate event to svg itself
            event.stopPropagation();
        }

        // attach click event to svg to return to default view
        chartSVG
            .attr("tabindex", 0)
            .on('click', () => {
                returnToDefaultView();
            })
            .on("keydown", (event, d) => {
                if (event.code == 'Enter' | event.code == 'Space') {
                    returnToDefaultView();
                }
            })

        function returnToDefaultView() {
            // Note that no country is selected
            countrySelected.value = false;
            // Reset y-axis label
            currentCountry.value = null;
            // arrange circles in DOM
            for (var i = 0; i < chartData.length; i++){
                let currentCircle = chartBounds.select("#circle-" + identifierAccessor(chartData[i]))
                let currentCircleMask = chartBounds.select("#circle-mask-" + identifierAccessor(chartData[i]))
                currentCircle.raise();
                currentCircleMask.raise();
            };
            chartBounds.selectAll(".circle")
                .attr('stroke', '#FFFFFF')
                .attr('stroke-width', '1px')
            chartBounds.selectAll(".circle-mask")
                .attr('fill-opacity', 0)
                .attr('stroke-opacity', 0)
            resetInfoBox()
        }

        const notationWidth = mobileView ? 160 : 160;
        chartSVG.append("text")
            .attr("aria-hidden", "true")
            .attr("class", "chart-notation")
            .attr("x", chartDimensions.boundedWidth - notationWidth / 2 )
            .attr("y", mobileView ? 2 : 0)
            .attr("dx", '0em')
            .attr("dy", '0em')
            .attr("dominant-baseline", "text-before-edge")
            .style("text-anchor", "middle")
            .attr("text-width", notationWidth)
            .text("More consumption of more vulnerable species")
            .call(d => wrap(d, {shift: false}));

        chartBounds.append("text")
            .attr("aria-hidden", "true")
            .attr("class", "chart-notation")
            .attr("x", mobileView ? 90 : 160)
            .attr("y", chartDimensions.boundedHeight * 0.99)
            .attr("dx", '0em')
            .attr("dy", mobileView ? '-1.75em' : '-2em')
            .attr("dominant-baseline", "text-after-edge")
            .style("text-anchor", "middle")
            .attr("text-width", notationWidth)
            .text("Less consumption of less vulnerable species")
            .call(d => wrap(d, {shift: false}));
    }

    function resetInfoBox() {
        activeCountry.value = defaultInfo;
    }

    function populateYAxisLabel(country, value) {
        currentCountry.value = country == 'United States of America' ? `the United States of America` : country;
        currentVal.value = value < 1 ? (value * 1000 < 10 ? `${Math.round(value * 10000) / 10} g` : `${Math.round(value * 1000 / 5) * 5} g`) : `${Math.round(value * 10) / 10} kg`;
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
            baseline = text.attr("dominant-baseline"),
            x = text.attr("x"),
            y = text.attr("y"),
            dy = parseFloat(text.attr("dy")),
            dx = parseFloat(text.attr("dx")),
            tspan = text.text(null).append("tspan").attr("y", y).attr("dy", dy + "em").attr("dominant-baseline", baseline);;

            while ((word = words.pop())) {
            line.push(word);
            tspan.text(line.join(" "));
                if (tspan.node().getComputedTextLength() > width) {
                line.pop();
                tspan.text(line.join(" "));
                line = [word];
                tspan = text.append("tspan").attr("x", x).attr("y", y).attr("dx", dx).attr("dy", ++lineNumber * lineHeight + dy + "em").attr("dominant-baseline", baseline).text(word);
                }
            }

            // https://stackoverflow.com/questions/60558291/wrapping-and-vertically-centering-text-using-d3-js
            if (lineNumber > 0  && shift) {
                const startDy = -(lineNumber * (lineHeight / 2));
                text
                    .selectAll("tspan")
                    .attr("dy", (d, i) => startDy + lineHeight * i + "em");
            }
        }
    )};
</script>

<style scoped lang="scss">
  #chart-container {
    width: 100%;
    margin: 0rem auto 2rem auto;
  }
  #grid-container {
    display: flex;
    flex-direction: column;
    margin: 3rem auto 4rem auto;
    width: 60vw;
    max-width: 1000px;
    @media only screen and (max-width: 600px) {
      width: 100%;
      flex-direction: column;
    }
  }
  .axis-title-container {
    margin-top: 4rem;
    margin-left: 0;
    height: 3.6rem;
    @media only screen and (max-width: 600px) {
      height: 5.4rem;
    }
  }
  .axis-title-container .axis-title {
    line-height: 1;
  }
</style>

<style lang="scss">
  /* css for elements added and classed w/ d3 */
  #wrapper {
    padding: 2px 0 0 2px;
  }
  #wrapper:focus {
    outline: none;
    border: none;
  }
  #wrapper:focus-visible {
    padding: -2px 0 0 -2px;
    border: none;
    outline: 2px solid #000000;
    outline-offset: 2px;
  }
  #y-axis line {
    stroke-width: 0.5px;
    stroke: var(--dark-grey);
    stroke-dasharray: 2 4;
  }
  .axis-text {
    font-size: 1.6rem;
    font-family: var(--default-font);
    user-select: none;
  }
  .axis-title {
    font-size: 1.8rem;
    font-family: var(--default-font);
    font-weight: 800;
    fill: var(--color-text);
    user-select: none;
  }

  p.axis-title {
    padding: 0;
  }
  .axis-subtitle {
    font-size: 1.8rem;
    font-family: var(--default-font);
    font-style: italic;
    font-weight: 500;
    fill: var(--color-text);
    user-select: none;
  }
  .chart-notation {
    font-size: 1.6rem;
    font-family: var(--default-font);
    font-style: italic;
    font-weight: 300;
    fill: var(--grey-12_1);
    user-select: none;
  }
  .circle-mask:focus-visible {
    stroke-opacity: 1;
    outline: none;
    border: none;
    stroke: #000000;
    stroke-width: 3px;
  }
  .circle-mask:focus {
    outline: none;
    border: none;
  }
  #consumable-tooltip {
    margin-left: -180px;
    width: 185px;
    z-index: 10;
    @media only screen and (max-width: 600px) {
        width: 160px;
        margin-left: -160px;
    }
    }
  #vulnerability-tooltip-span {
    z-index: 9;
  }
  #vulnerability-tooltip {
    margin-left: -235px;
    width: 235px;
    @media only screen and (max-width: 600px) {
        width: 210px;
        margin-left: -210px;
    }
  }
  #vulnerability-tooltip a {
    color: white;
    text-decoration: underline;
    font-weight: 400;
  }
</style>