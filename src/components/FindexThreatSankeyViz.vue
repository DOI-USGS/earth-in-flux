<template>
    <!---VizSection-->
    <VizSection
        :figures="true"
        :fig-caption="false"
    >
        <!-- HEADING -->
        <template #heading>
            <h2>
            </h2>
        </template>
        <!-- FIGURES -->
        <template #aboveExplanation>
                <p v-html="text.paragraph1" />
        </template>
        <template #figures>
            <div id="threat-container" ref="chart"></div>
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
    import * as d3sankey from 'd3-sankey';
    import VizSection from '@/components/VizSection.vue';

    // define props
    defineProps({
        text: { type: Object }
    })

    // Global variables 
    const publicPath = import.meta.env.BASE_URL;
    const mobileView = isMobile;
    const data = ref();
    const dataFile = 'findex_total_weighted_threats.csv'
    const chart = ref(null);
    let chartDimensions;
    let chartBounds;
    let nodeGroup;
    let linkGroup;
    let textGroup;
    const nodeWidth = 4;
    const labelBuffer = 10;

    // Colors for threat categories, Needs to be updated with CSS for text legend
    const categoryColors = {
        'Climate and weather': '#899bb7',
        'Fishing pressure': '#c29fcd',
        'Habitat': '#729C9D', 
        'Invasive species': '#d38884',
        'Pollution': '#dab589'
    }; 

    onMounted(async () => {
        try {
            await loadDatasets();
            
            if (data.value.length > 0) {
                initSankey({
                    width: chart.value.offsetWidth,
                    height: window.innerHeight * 0.8,
                    margin: 10,
                    marginLeft: mobileView ? 80: 150,
                    marginRight: mobileView ? 135: 250,
                    marginTop: mobileView ? 60 : 50,
                    containerId: 'threat-container'
                });
                createSankey({
                    dataset: data.value,
                    containerId: 'threat-container'
                });
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
    };

    async function loadData(fileName) {
        try {
            const data = await d3.csv(publicPath + fileName, d => {
                return d;
            });
            return data;
        } catch (error) {
            console.error(`Error loading data from ${fileName}`, error);
            return [];
        }
    };

    function initSankey({
        width,
        height,
        margin,
        marginTop = margin, // top margin, in pixels
        marginBottom = margin, // left margin, in pixels
        marginLeft = margin, // left margin, in pixels
        marginRight = margin, // right margin, in pixels
        containerId
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


        // draw svg canvas for sankey
        const svg = d3.select('#' + containerId)
            .append('svg')
            .attr('class', 'sankeySVG')
            .attr('viewBox', `0 0 ${chartDimensions.width} ${chartDimensions.height}`)
            .style('width', "100%")
            .style('height', "100%");

        // add group for bar chart bounds, translating by chart margins
        chartBounds = svg.append('g')
            .attr('id', 'wrapper')
            .style("transform", `translate(${
                chartDimensions.margin.left
            }px, ${
                chartDimensions.margin.top
            }px)`)

        // Add group to chart bounds to hold all sankey path groups
        nodeGroup = chartBounds.append('g')
            .attr('id', 'node_group')
        
        linkGroup = chartBounds.append('g')
            .attr('id', 'link_group')

        textGroup = chartBounds.append('g')
            .attr('id', 'text_group')

        // add titles
        const subtitle = mobileView ? "High to low" : "Ranked high to low"
        const leftTitle = svg.append("text")
            .attr("class", "axis-title")
            .attr("x", chartDimensions.margin.left - labelBuffer + nodeWidth) // match spacing between sankey and labels
            .attr("y", 0)
            .attr("dx", "0em")
            .attr("dy", "0em")
            .attr("dominant-baseline", "hanging")
            .attr("text-width", chartDimensions.margin.left)
            .style("text-anchor", "end")
            .text("Threat Categories")
            .call(d => mobileView ? wrap(d, {shift: false}) : d)

        const leftTitleLength = leftTitle.node().getComputedTextLength()

        svg.append("text")
            .attr("class", "axis-text axis-value axis-notation")
            .attr("x", chartDimensions.margin.left - labelBuffer + nodeWidth) // match spacing between sankey and labels
            .attr("y", 0)
            .attr("dx", "0em")
            .attr("dy", leftTitleLength > chartDimensions.margin.left ? "2.8em" : "1.5em")
            .attr("dominant-baseline", "hanging")
            .attr("text-width", chartDimensions.margin.left)
            .attr("text-anchor", "end")
            .text(subtitle)
            .call(d => wrap(d, {shift: false}))

        const rightTitle = svg.append("text")
            .attr("class", "axis-title")
            .attr("x", chartDimensions.width - chartDimensions.margin.right + labelBuffer - nodeWidth) // match spacing between sankey and labels
            .attr("y", 0)
            .attr("dx", "0em")
            .attr("dy", "0em")
            .attr("dominant-baseline", "hanging")
            .attr("text-width", chartDimensions.margin.right)
            .style("text-anchor", "start")
            .text("Threats")

        const rightTitleLength = rightTitle.node().getComputedTextLength()

        svg.append("text")
            .attr("class", "axis-text axis-value axis-notation")
            .attr("x", chartDimensions.width - chartDimensions.margin.right + labelBuffer - nodeWidth) // match spacing between sankey and labels
            .attr("y", 0)
            .attr("dx", "0em")
            .attr("dy", rightTitleLength > chartDimensions.margin.right ? "2.8em" : "1.5em")
            .attr("dominant-baseline", "hanging")
            .attr("text-width", chartDimensions.margin.right)
            .attr("text-anchor", "start")
            .text(subtitle)
            .call(d => wrap(d, {shift: false}))
    };

    function createSankey({
        dataset
    }) {

        // get unique categories and parameters
        const categoryGroups = [... new Set(dataset.map(d => d.ThreatCategory))];
    
        // initialize sankey
        const sankey = d3sankey.sankey()
            .nodeSort((a,b) => d3.descending(a.value, b.value))
            .linkSort(null)
            .nodeWidth(nodeWidth)
            .nodePadding(mobileView ? 15 : 11)
            .extent([[0, 0], [chartDimensions.boundedWidth, chartDimensions.boundedHeight]])

        // Set up color scale 
        const colorScale = d3.scaleOrdinal()
            .domain(categoryGroups)
            .range(categoryGroups.map(item => categoryColors[item]));
        
        // set up the nodes and links
        var nodesLinks = graphNodes({
            data: dataset
        });

        const {nodes, links} = sankey({
            nodes: nodesLinks.nodes.map(d => Object.create(d)),
            links: nodesLinks.links.map(d => Object.create(d))
        });

        // Set up transition.
        const dur = 1000;

        // Update nodes for sankey, assigning data
        nodeGroup.selectAll('g')
            .data(nodes)
            .join(
                enter => enter
                    .append('rect')
                        .attr("x", d => d.x0)
                        .attr("y", d => d.y0)
                        .attr("height", d => d.y1 - d.y0)
                        .attr("width", d => d.x1 - d.x0)
                    .append("title")
                        .text(d => `${d.name}\n${d.value.toLocaleString()}`),

                null, // no update function

                exit => {
                    exit
                    .transition()
                    .duration(dur / 2)
                    .style("fill-opacity", 0)
                    .remove();
            });

        // Update links for sankey, assigning data
        linkGroup.selectAll('g')
            .data(links)
            .join(
                enter => {
                    enter 
                        .append("path")
                            .attr("d", d3sankey.sankeyLinkHorizontal())
                            .attr("stroke", d => colorScale(d.names[0]))
                            .attr("stroke-width", d => d.width)
                            .style("mix-blend-mode", "multiply")
                            .style('fill', "none")
                        .append("title")
                            .text(d => `${d.names.join(" → ")}\n${d.value.toLocaleString()}`)
                },

                null,

                exit => {
                    exit
                        .transition()
                        .duration(dur / 2)
                        .style("fill-opacity", 0)
                        .style("stroke-width", 0)
                        .style("color-opacity", 0)
                        .remove();
                }
            );


        // Update text for sankey, assigning data from nodes
        textGroup.selectAll('g')
            .data(nodes)
            .join(
                enter => {
                    enter
                        .append("text")
                            .attr("class", d => d.x0 < chartDimensions.boundedWidth / 2 ? "axis-text left" : "axis-text right")
                            .attr("x", d => d.x0 < chartDimensions.boundedWidth / 2 ? d.x1 : d.x0) //checks for right-most labels
                            .attr("y", d => (d.y1 + d.y0) / 2)
                            .attr("dy", "0em")
                            .attr("dx", d => d.x0 < chartDimensions.boundedWidth / 2 ? -labelBuffer : labelBuffer)
                            .attr("dominant-baseline", "central")
                            .attr("text-anchor", d => d.x0 < chartDimensions.boundedWidth / 2 ? "end" : "start") //checks for right-most labels
                            .attr("text-width", d => d.x0 < chartDimensions.boundedWidth / 2 ? chartDimensions.margin.left : chartDimensions.margin.right)
                            .text(d => d.name)
                            .call(d => mobileView ? wrap(d, {shift: true}) : d)
                        // .append("tspan")
                        //     .attr("fill-opacity", 0.7)
                        //     .text(d => ` ${d.value.toLocaleString()}`)
                        //     .style("font", "14px sans-serif")
                },
                null,
                exit => {
                    exit
                        .transition()
                        .duration(dur / 2)
                        .style("fill-opacity", 0)
                        .style("stroke-width", 0)
                        .style("color-opacity", 0)
                        .remove();
                }
            );
    };

    // set up the nodes and links
    function graphNodes({data}){ //https://observablehq.com/@d3/parallel-sets?collection=@d3/d3-sankey
        let keys = data.columns.slice(0, 2); // which columns for nodes

        let index = -1;
        let nodes = [];
        let nodeByKey = new d3.InternMap([], JSON.stringify);
        let indexByKey = new d3.InternMap([], JSON.stringify);
        let links = [];

        // creates nodes for each column (keys)
        for (const k of keys) {
            for (const d of data) {
            const key = [k, d[k]];
            if (nodeByKey.has(key)) continue;
                const node = {name: d[k]};
                nodes.push(node);
                nodeByKey.set(key, node);
                indexByKey.set(key, ++index);
            }
        }

        // creates links between nodes
        for (let i = 1; i < keys.length; ++i) {
            const a = keys[i - 1];
            const b = keys[i];
            const prefix = keys.slice(0, i + 1);
            const linkByKey = new d3.InternMap([], JSON.stringify);
            for (const d of data) {
            const names = prefix.map(k => d[k]);
            const value = d.TotalWeightedThreatMetric; 
            let link = linkByKey.get(names);
            if (link) { link.value += value; continue; }
            link = {
                source: indexByKey.get([a, d[a]]),
                target: indexByKey.get([b, d[b]]),
                names,
                value
            };
            links.push(link);
            linkByKey.set(names, link);
            }
        }
        return {nodes, links};
    };

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
</script>

<style lang="scss">
    #threat-container {
        max-width: min(1000px, 60vw);
        margin: 5rem auto 0 auto;
        @media screen and (max-width: 600px) {
            max-width: 100%;
        }
    }
    .axis-text {
        font-size: 1.6rem;
        font-family: var(--default-font);
        user-select: none;
        @media screen and (max-width: 600px) {
            font-size: 1.4rem;
        }
    }
    .axis-title {
        font-size: 1.8rem;
        font-family: var(--default-font);
        font-weight: 900;
        fill: var(--color-text);
        user-select: none;
        @media screen and (max-width: 600px) {
            font-size: 1.6rem;
        }
    }
    .axis-notation {
        font-style: italic;
    }
</style>