<template>
    <!---VizSection-->
    <VizSection
        id="cross-section"
        :figures="true"
        :fig-caption="false"
    >
        <!-- HEADING -->
        <template #heading>
            <h2>
                {{ text.heading }}
            </h2>
        </template>
        <!-- FIGURES -->
        <template #aboveExplanation>
            <p v-html="text.paragraph1" />
        </template>
        <template #figures>
            <div id="aerosols-grid-container" />
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
    import { onMounted } from "vue";
    import * as d3 from 'd3';
    import VizSection from '@/components/VizSection.vue';

    // define props
    defineProps({
        text: { type: Object }
    })

    // Declare behavior on mounted
    // functions called here
    onMounted(async () => {
        try {
            // Use external svg from s3
            d3.xml("https://labs.waterdata.usgs.gov/visualizations/svgs/aerosols.svg").then(function(xml) {
                // add svg content to DOM
                const svgGrid = document.getElementById("aerosols-grid-container")
                svgGrid.appendChild(xml.documentElement);
                
                // add id to svg
                d3.select("#aerosols-grid-container").select("svg")
                    .attr("id", "aerosols-svg")

                // add interactivity
                addInteractions();
            });
        } catch (error) {
            console.error('Error during component mounting', error);
        }        
    });


    function draw_trajectories(fire){
        for(let trajectory_num=1;trajectory_num<=24;trajectory_num++){
            d3.select("#trajectory-" + fire + "-" + trajectory_num).selectAll("path")
                .style("stroke-opacity", 0.4);
            var trajectory_line = d3.select("#trajectory-" + fire + "-" + trajectory_num).selectAll("path");
            var totalLength = trajectory_line.node().getTotalLength();
            // Animate the line
            trajectory_line.style("transition", "none")
                .style("stroke-dasharray", totalLength + " " + totalLength)
                .style("stroke-dashoffset", totalLength)
                .transition()
                .delay(200*(trajectory_num-1))
                .duration(2000)
                .ease(d3.easeLinear)
                .style("stroke-dashoffset", 0);

            d3.select('#wildfire-label-'+fire).select("text")
                .style("opacity", 1.0);
        }
    }

    function remove_trajectories(fire){
        for(let trajectory_num=1;trajectory_num<=24;trajectory_num++){
            d3.select("#trajectory-" + fire + "-" + trajectory_num).selectAll("path")
                .style("stroke-opacity", 0.0);
            d3.select('#wildfire-label-'+fire).select("text")
                .style("opacity", 0.0);
        }
    }

    function highlight_core(coreid, corealpha, corecolor){
        console.log(coreid)
        d3.select("#"+coreid).selectAll("path")
            .style("stroke-opacity", corealpha)
            .style("stroke", corecolor);
    }

    function mouseover(event,coreids) {
        if (event.currentTarget.id.startsWith("source-")){
            let fire = event.currentTarget.id.slice(7);
            draw_trajectories(fire);
            highlight_core(coreids[fire], 1.0, "#d62728");
        }
    }

    function mouseout(event,coreids) {
        if (event.currentTarget.id.startsWith("source-")){
            let fire = event.currentTarget.id.slice(7);
            remove_trajectories(fire);
            highlight_core(coreids[fire], 0.2, "#000000");
        }
    }

    function addInteractions() {
        // set viewbox for svg with loss function chart
        const aerosolsSVG = d3.select("#aerosols-svg");

        const coreids = ["id-2017-Core-1","id-2016-Core-3"]
        // Add interaction to loss function chart
        aerosolsSVG.selectAll("g")
            .on("mouseover", (event) => mouseover(event,coreids))
            .on("mouseout", (event) => mouseout(event,coreids));
    }
</script>

<style scoped lang="scss">
    #aerosols-grid-container {
        display: grid;
        width: 100%;
        max-width: 1200px;
        margin: 0 auto 0 auto;
        grid-template-areas:
            "chart";
    }
</style>
<style lang="scss">
/* css for elements added/classed w/ d3 */
    #aerosols-svg {
        grid-area: chart;
        place-self: center;
        max-height: 100%;
        max-width: 100%;
    }
</style>