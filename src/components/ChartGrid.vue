<template>
    <div id="chartGrid" class="padded">
        <a v-for="(item, index) in chartContent" :key="index" :href="item.route">
            <div class="chart">
                <img :src="item.img_src" alt=""/>
            </div>
        </a>
    </div>
</template>

<script>
import ChartGrid from '@/assets/content/ChartGrid.js';

export default {
  data() {
    return {
      chartContent: ChartGrid.chartGridItems
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