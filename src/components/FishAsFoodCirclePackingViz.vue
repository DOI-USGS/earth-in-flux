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
        <<template #aboveExplanation>
        <p v-html="text.paragraph1" />
        <FamilyInfoBox :activeFamily="activeFamily" />
        </template>

        <template #figures>
            <div class="chart-container single" ref="chart">
            </div>
            
        </template>
        <!-- FIGURE CAPTION -->
        <template #figureCaption>
        </template>
        <!-- EXPLANATION -->
        <template #belowExplanation>
            <p v-html="text.paragraph2" />
        </template>
    </VizSection>
</template>

<script setup>
    import { onMounted, ref } from "vue";
    import * as d3 from 'd3';
    import VizSection from '@/components/VizSection.vue';
    import FamilyInfoBox from '@/components/FamilyInfoBox.vue';
    

    // define props
    defineProps({
        text: { type: Object }
    })

    // global variables
    const publicPath = import.meta.env.BASE_URL;
    const chart = ref(null);
    const bodyCSS = window.getComputedStyle(document.body);
    const lightBlue = bodyCSS.getPropertyValue('--light-blue');
    const darkBlue = bodyCSS.getPropertyValue('--dark-blue');

    const familyInfo = {
        "Cyprinidae": {
            image: "https://labs.waterdata.usgs.gov/visualizations/images/common-carp.jpeg",
            text: "Cyprinidae is a family of freshwater fish commonly called the carp or minnow family, including the carps, the true minnows, and their relatives the barbs and barbels, among others. Cyprinidae is the largest and most diverse fish family, and the largest vertebrate animal family overall, with about 1,780 species divided into 166 valid genera..... ",
            caption: "Common carp <em>(Cyprinus carpio)</em>"
        },
        "Salmonidae": {
            image: "https://labs.waterdata.usgs.gov/visualizations/images/pink-salmon.jpeg",
            text: "Salmonidae is a family of ray-finned fish that constitutes the only currently extant family in the order Salmoniformes, consisting of 11 extant genera and over 200 species collectively known as salmonids or salmonoids. The family includes salmon (both Atlantic and Pacific species), trout (both ocean-going and landlocked), char, graylings, freshwater whitefishes, taimens and lenoks, all coldwater mid-level predatory fish that inhabit the subarctic and cool temperate waters of the Northern Hemisphere...",
            caption: "Pink salmon <em>(Oncorhynchus gorbuscha)</em>"
        },
        "Centrarchidae": {
            image: "https://labs.waterdata.usgs.gov/visualizations/images/black-crappie.jpeg",
            text: "Centrarchidae, better known as sunfishes, is a family of freshwater ray-finned fish belonging to the order Centrarchiformes, native only to North America. The centrarchid family comprises 38 identified species, 34 of which are extant. It includes many popular game fishes familiar to North American anglers, such as the rock bass, largemouth bass, bluegill, pumpkinseed, green sunfish and crappies... ",
            caption: "Black crappie <em>(Pomoxis nigromaculatus)</em>"
        }
        };

    const defaultFamily = {
            image: "https://labs.waterdata.usgs.gov/visualizations/images/default-fish.jpeg", 
            text: "Click on a fish family in the chart to learn more about its importance and characteristics.",
            caption: "Brown trout <em>(Salmo trutta)</em>"
        };

const activeFamily = ref(defaultFamily); // start with placeholder
    
    // Declare behavior on mounted
    // functions called here
    onMounted(async () => {
        // Load the json data
        const data = await d3.json(publicPath + 'total_price.json');

        // build chart
        buildChart(data);
    });

    // https://observablehq.com/@d3/zoomable-circle-packing
    function buildChart(data) {

        // Specify the chart’s dimensions.
        const width = 928;
        const height = width;

        // Create the color scale.
        const color = d3.scaleLinear()
            .domain([0, 5])
            .range([lightBlue, darkBlue])
            .interpolate(d3.interpolateHcl);

        // Compute the layout.
        const pack = data => d3.pack()
            .size([width, height])
            .padding(3)
        (d3.hierarchy(data)
            .sum(d => d.value)
            .sort((a, b) => b.value - a.value));
        const root = pack(data);

        // Create the SVG container.
        const svg = d3.select(chart.value)
            .append("svg")
            .attr("id", "circle-pack-svg")
            .attr("viewBox", `-${width / 2} -${height / 2} ${width} ${height}`)
            .attr("width", width)
            .attr("height", height)

        // Append the nodes.
        const node = svg.append("g")
            .selectAll("circle")
            .data(root.descendants().slice(1))
            .join("circle")
                .attr("fill", d => d.children ? color(d.depth) : "white")
                .attr("pointer-events", d => !d.children ? "none" : null)
                .on("mouseover", function() { d3.select(this).attr("stroke", "#000"); })
                .on("mouseout", function() { d3.select(this).attr("stroke", null); })
                .on("click", (event, d) => {
                    if (focus !== d) {
                        zoom(event, d);
                        event.stopPropagation();
                    }
                    });


        // Append the text labels.
        const label = svg.append("g")
            .attr("class", "circle-packing-text")
            .attr("pointer-events", "none")
            .attr("text-anchor", "middle")
            .selectAll("text")
            .data(root.descendants())
            .join("text")
                .attr("class", "fish-title")
                    .style("fill-opacity", d => d.parent === root ? 1 : 0)
                    .style("stroke-opacity", d => d.parent === root ? 1 : 0)
                    .style("display", d => d.parent === root ? "inline" : "none")
                    .style("font-size", "1.2rem")
                    .style("stroke", "white")          // white outline
                    .style("stroke-width", "2px")      
                    .style("paint-order", "stroke")    
                    .style("stroke-linejoin", "round") 
                    .style("fill", "black")            
                    .text(d => [d.data.name, d3.format("$.1s")(d.value).replace("G","B")].join("\n"));


        // Create the zoom behavior and zoom immediately in to the initial focus node.
        svg.on("click", (event, d) => {
            if (focus !== d) {
                zoom(event, d);
                event.stopPropagation();
                }
            });

        //  zoom-out on background click
        d3.select(chart.value).select("svg").on("click", (event) => {
            activeFamily.value = defaultFamily;
            zoom(event, root);
        });

        let focus = root;
        let view;
        zoomTo([focus.x, focus.y, focus.r * 2.2]);

        function zoomTo(v) {
            const k = width / v[2];

            view = v;

            label.attr("transform", d => `translate(${(d.x - v[0]) * k},${(d.y - v[1]) * k})`);
            node.attr("transform", d => `translate(${(d.x - v[0]) * k},${(d.y - v[1]) * k})`);
            node.attr("r", d => d.r * k);
        }

        function zoom(event, d) {
                focus = d;

                const transition = svg.transition()
                    .duration(event.altKey ? 7500 : 750)
                    .tween("zoom", () => {
                    const i = d3.interpolateZoom(view, [focus.x, focus.y, focus.r * 2.2]);
                    return t => zoomTo(i(t));
                    });

                label
                    .filter(function (d) {
                    return d.parent === focus || this.style.display === "inline";
                    })
                    .transition(transition)
                    .style("fill-opacity", d => d.parent === focus ? 1 : 0)
                    .style("stroke-opacity", d => d.parent === focus ? 1 : 0)
                    .on("start", function (d) {
                    if (d.parent === focus) this.style.display = "inline";
                    })
                    .on("end", function (d) {
                    if (d.parent !== focus) this.style.display = "none";
                    });

                    let current = d;
                    let foundFamily = null;

                    while (current) {
                    const name = current.data?.name;
                    if (name && familyInfo[name]) {
                        foundFamily = name;
                        break;
                    }
                    current = current.parent;
                    }

                    if (foundFamily) {
                    activeFamily.value = familyInfo[foundFamily];
                    } else {
                    activeFamily.value = defaultFamily;
                    }

            }

        return svg.node();
    }

</script>

<style lang="scss">
/* css for elements added/classed w/ d3 */
    #circle-pack-svg {
        max-width: 100%;
        height: auto;
    }
    .circle-packing-text {
        font-size: 1.2rem;
        @media screen and (max-width: 600px) {
            font-size: 2rem;
        }
    }
</style>
