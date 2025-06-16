<template>
    <!---VizSection-->
    <VizSection
        :figures="true"
        :fig-caption="false"
    >
        <template #figures>
          <div
            id="grid-container"
          >
            <CountryInfoBox :activeCountry="activeCountry"  />
            <div
              id="sub-grid-container"
            >
                <div
                    id="continent-map-image-container"
                >
                    <img 
                    id="continent-map-image"
                    :src="getImageURL(text.image)"
                    :alt="text.imageAlt"
                    >
                </div>
                <div id="chart-container" class="maxWidth" ref="chart" />           
            </div>
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
    import { isMobile } from 'mobile-device-detect';
    import * as d3 from 'd3';
    import VizSection from '@/components/VizSection.vue';
    import CountryInfoBox from '@/components/CountryInfoBox.vue';

    // define props
    const props = defineProps({
        text: { type: Object }
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
    const chartTitle = 'Title of chart';
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
        Africa: "#648E8E",
        Oceania: "#845c93",
        "North America": "#b24d4b",
        Europe: "#a27846",
        Asia: "#1D3867",
        "South America": "#899bb7"
    }
    const colorLowerBound = "#6F4E37" // lower end
    const colorUpperBound =  "#5CB270" // upper end
    let colorScale;

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
                    marginBottom: 60,
                    marginLeft:5,
                    marginTop: mobileView ? 80 : 30
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
            console.log(`${file} data in`);
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
        customFormat = false,
        tickValues = null,
        textAngle = 0,
        keepDomain = true,
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
        customFormat = false,
        tickValues = null,
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
            customFormat: customFormat,
            tickValues: tickValues,
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
        customFormat = false,
        tickValues = null,
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
            customFormat: customFormat,
            tickValues: tickValues,
            keepDomain: keepDomain,
        })
    }

    function initRScale() {
        rScale = d3.scaleSqrt()
            .range([Math.min(chartDimensions.boundedWidth,chartDimensions.boundedHeight)*rPropMin,Math.min(chartDimensions.boundedWidth,chartDimensions.boundedHeight)*rPropMax]);
    }

    function initColorScale() {
        colorScale = d3.scaleLinear()            
            .range([colorLowerBound, colorUpperBound]);
    }

    function drawChart(data) {
        //////////////////////////////
        /////    PROCESS DATA    /////
        //////////////////////////////
        const chartData = data //.filter(d => d.continent === continent); // May filter later

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
        const xInnerDomainRange = d3.max(chartData, xAccessor) - d3.min(chartData, xAccessor);
        let xDomain_min = d3.min(chartData, xAccessor) - xInnerDomainRange * rPropMax / (1.0 - 2.0 * rPropMax); //ensures that the buffer is the max radius away from the smallest data point
        xDomain_min = mobileView ? xDomain_min * 0.95 : xDomain_min;
        const xDomain_max = d3.max(chartData, xAccessor) + xInnerDomainRange * rPropMax / (1.0 - 2.0 * rPropMax); //ensures that the buffer is the max radius away from the larger data point        
        xScale
            .domain([xDomain_min, xDomain_max]);
        drawXAxis({axisTitle: 'Climate vulnerability', tickValues: []})
        //console.log(r_prop_max * (xDomain_max - xDomain_max), xInnerDomainRange * r_prop_max / (1.0 - 2.0 * r_prop_max))
        // Add arrow
        const arrow_dim = 12;
        chartSVG.append("defs").append("marker")
            .attr("id", "arrowhead-right") // Assign a unique ID
            .attr("viewBox", "0 -5 10 10") // Define the coordinate system for the marker
            .attr("refX", 10) // Reference point for the end of the arrow
            .attr("refY", 0)
            .attr("markerWidth", arrow_dim)
            .attr("markerHeight", arrow_dim)
            .append("path")
            .attr("d", "M 0,-5 L 10,0 L 0,5"); // Define the arrowhead shape
        chartSVG.append("defs").append("marker")
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
            .attr("y", chartDimensions.margin.bottom / 2)
            .attr("dominant-baseline", "text-after-edge")
            .style("text-anchor", "start")
            .text("lower​")

        d3.select("#x-axis").append("text")
            .attr("class", "axis-subtitle")
            .attr("x", chartDimensions.boundedWidth)
            .attr("y", chartDimensions.margin.bottom / 2)
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
              suffix = "g"              
              return d3.format(".1s")(d).replace("m","") + " " + suffix;
            } else {
              suffix = "kg"              
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

        chartSVG.append("text")
            .attr("class", "axis-title")
            .attr("x", chartDimensions.margin.left)
            .attr("y", 0)
            .attr("dx", mobileView ? "0em" : "-0.2em")
            .attr("dy", "0em")
            .attr("dominant-baseline", "hanging")
            .attr("text-width", chartDimensions.boundedWidth)
            .style("text-anchor", "start")
            .text("Per fisher consumption of recreationally-harvested inland fish is…​")
            .call(d => mobileView ? wrap(d, {shift: false}) : d)

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
        initColorScale()
        colorScale
            .domain(d3.extent(chartData, colorAccessor));

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
                    .attr('fill', d => colors[colorAccessor(d)])
                    .attr('stroke', '#FFFFFF')
                    // .on("mouseover", (event, d) => {
                    //     activeCountry.value = {
                    //         name: d.admin
                    //     }
                    // })
                    .on("click", (event, d) => {
                        console.log(d.admin)
                        activeCountry.value = {
                            name: d.admin,
                            population: d3.format(',')(d.population),
                            participation_rate: Math.round(d.participation_rate),
                            n_fishers: d3.format(',')(d.n_fishers),
                            consum_harv_kg: d3.format(',')(d.total_consumable_harv_kg),
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
                                    percent: Math.round(d.unkown*10)/10
                                }
                            ]
                        }
                    })
                    // .on("mouseout", resetInfoBox)
        d3.select(chart.value)
            .select('svg')
            .on('click', () => {
              resetInfoBox
            })
    }

    function resetInfoBox() {
        console.log('reset')
        console.log(defaultInfo)
        activeCountry.value = defaultInfo;
        console.log(activeCountry.value)
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
            
            console.log(text.attr("dy"))

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

    function getImageURL(filename) {
        return new URL(`../assets/images/${filename}`, import.meta.url).href
    }
</script>

<style scoped lang="scss">
  #chart-container {
    width: 60vw;
    @media only screen and (max-width: 600px) {
      width: 100%;
    }
  }
  #grid-container {
    display: flex;
    flex-direction: column;
    margin: 3rem auto 4rem auto;
    @media screen and (max-width: 600px) {
        flex-direction: column;
    }
  }
  #sub-grid-container {
    width: 60vw;
    display: flex;
    flex-direction: row;
    position: relative;
    margin: 6rem auto 2rem auto;
    @media screen and (max-width: 600px) {
        width: 100%;
        flex-direction: column;
    }
  }
  #continent-map-image-container {
    position: absolute;
    right: 10px;
    @media screen and (max-width: 600px) {
        position: relative;
    }
  }
  #continent-map-image {
    width: 300px;
    @media screen and (max-width: 600px) {
        margin-top: 1rem;
    }
  }
</style>

<style lang="scss">
/* css for elements added and classed w/ d3 */
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
    font-weight: 900;
    fill: var(--color-text);
    user-select: none;
  }
  .axis-subtitle {
    font-size: 1.8rem;
    font-family: var(--default-font);
    font-style: italic;
    font-weight: 500;
    fill: var(--color-text);
    user-select: none;
  }
</style>