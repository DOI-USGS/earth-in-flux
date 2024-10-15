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
    import * as d3 from 'd3';
    import * as d3sankey from 'd3-sankey';
    import VizSection from '@/components/VizSection.vue';
import { isYandex } from "mobile-device-detect";

    // define props
    defineProps({
        text: { type: Object }
    })

    // Global variables 
    const publicPath = import.meta.env.BASE_URL;
    const data = ref();
    const dataFile = 'findex_total_weighted_threats.csv'
    const chart = ref(null);
    const chartTitle = 'Title of chart';
    let chartDimensions;
    let chartBounds;
    let nodeGroup;
    let linkGroup;
    let textGroup;

    // Colors for threat categories, Needs to be updated with CSS for text legend
    const categoryColors = {
        'Climate and weather': '#EECEB9',
        'Exploitation':  '#939185',
        'Habitat':  '#C8ACD6', 
        'Invasive species':  '#80909D',
        'Pollution': '#E8E8E3'
    }; 

    onMounted(async () => {
        try {
            await loadDatasets();
            
            if (data.value.length > 0) {
                initSankey({
                    width: chart.value.offsetWidth,
                    height: window.innerHeight * 0.8,
                    margin: 10,
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
                // d.TotalWeightedThreatMetric = + d.TotalWeightedThreatMetric,
                // d.TotalWeightedThreatMetricx1000 = + d.TotalWeightedThreatMetricx1000
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
    };

    function createSankey({
        dataset,
        containerId
    }) {

        // get unique categories and parameters
        const categoryGroups = [... new Set(dataset.map(d => d.ThreatCategory))];
        // const categoryGroups = d3.union(d3.map(dataset, d => d.Category));
        const habitatGroups = d3.union(d3.map(dataset, d => d.Habitat));
    
        // initialize sankey
        const sankey = d3sankey.sankey()
            .nodeSort(null)
            .linkSort(null)
            .nodeWidth(4)
            .nodePadding(11)
            .extent([[150, 5], [chartDimensions.boundedWidth - 300, chartDimensions.boundedHeight - 0]])

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
        const t = d3.transition().duration(dur);

        // Update nodes for sankey, assigning data
        const sankeyNodesGroups = nodeGroup.selectAll('g')
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
        const sankeyLinksGroups = linkGroup.selectAll('g')
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
        const sankeyTextGroups = textGroup.selectAll('g')
            .data(nodes)
            .join(
                enter => {
                    enter
                        .append("text")
                            .attr("x", d => d.x0 < chartDimensions.boundedWidth / 2 ? d.x1 -10 : d.x0 + 10) //checks for right-most labels
                            .attr("y", d => (d.y1 + d.y0) / 2)
                            .attr("dy", "0.35em")
                            .attr("text-anchor", d => d.x0 < chartDimensions.boundedWidth / 2 ? "end" : "start") //checks for right-most labels
                            .text(d => d.name)
                            .style("font", "14px sans-serif")
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
                // Doing some custom index setting to sort the left side of the sankey
                if (d[k] == 'Invasive species') {
                    ++index // still need to advance index
                    indexByKey.set(key, 3); // Would otherwise be 1
                } else if (d[k] == 'Pollution') {
                    ++index
                    indexByKey.set(key, 1); // Would otherwise be 2
                } else if (d[k] == 'Climate and weather') {
                    ++index
                    indexByKey.set(key, 2); // Would otherwise be 3
                } else {
                    indexByKey.set(key, ++index);
                }
                // Use below if dropping custom sorting
                // indexByKey.set(key, ++index);
            }
        }
        
        // With custom indices, need to re-order nodes
        // Here, taking 'Invasive species' out of slot 1 and dropping in slot 3
        nodes.splice(3, 0, nodes.splice(1, 1)[0]);

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
</script>

<style>
</style>