<template>
    <div id="chartGrid" class="padded">
        <!-- When we get actual charts, we can pull all the info and urls in from a js file and loop through
        like in the VizLab home site -->
            <a href="/"> <!-- Populate hrefs with relative routes -->
                <div class="chart">
                    <img :src="sampleChart" alt=""/>
                </div>
            </a>
            <a href="/">
                <div class="chart">
                    <img :src="sampleChart" alt=""/>
                </div>
            </a>
            <a href="/">
                <div class="chart">
                    <img :src="sampleChart" alt=""/>
                </div>
            </a>
            <a href="/">
                <div class="chart">
                    <img :src="sampleChart" alt=""/>
                </div>
            </a>
            <a href="/">
                <div class="chart">
                    <img :src="sampleChart" alt=""/>
                </div>
            </a>
            <a href="/">
                <div class="chart">
                    <img :src="sampleChart" alt=""/>
                </div>
            </a>
            <a href="/">
                <div class="chart">
                    <img :src="sampleChart" alt=""/>
                </div>
            </a>
            <a href="/">
                <div class="chart">
                    <img :src="sampleChart" alt=""/>
                </div>
            </a>
    </div>
</template>

<script>
import sampleChart from '@/assets/images/sampleChart.PNG';

export default {
    data() {
        return {
            sampleChart
        };
    }
};
</script>

<style>

#chartGrid{
	margin: 50px auto;
	display: flex;
	justify-content: center;
	gap: 50px;
	flex-wrap: wrap;
	gap: 50px;
	width: 100%;
	max-width: 60%;

	.chart{
		display: flex;
		flex-direction: column-reverse;
		justify-content: flex-start;
		position: relative;
		width: 280px;
		height: 100%;
	}
    .chart:hover {
        transform: scale(1.1);
    }
}
</style>