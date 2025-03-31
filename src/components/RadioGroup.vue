<template>
  <div class="radio_wrap">
    <label
      v-for="option in options"
      :key="option.value"
      class="radio-label"
      :class="{ selected: modelValue === option.value }"
      :style="{
        backgroundColor: modelValue === option.value
          ? `${option.color}20`  // 20 is ~12% opacity in hex
          : 'transparent'
      }"
    >
      <input
        type="radio"
        class="radio-input"
        :value="option.value"
        :checked="modelValue === option.value"
        @change="$emit('update:modelValue', option.value)"
        name="radio-group"
      />
      <span
        class="radio-button"
        :style="{
          borderColor: modelValue === option.value ? option.color || activeColor : inactiveColor,
          backgroundColor: modelValue === option.value ? option.color || activeColor : 'transparent'
        }"
      ></span>

      <span
        class="radio-text"
        :class="{ ractive: modelValue === option.value }"
      >
        {{ option.label }}
      </span>
    </label>
  </div>
</template>

<script setup>
defineProps({
  modelValue: [String, Number], // v-model binding for selected value

  // array of radio options: [{ label: 'Option 1', value: 'opt1' }, ...]
  options: {
    type: Array,
    required: true
  },

  // optional colors
  activeColor: {
    type: String,
    default: 'var(--black-soft)'
  },
  inactiveColor: {
    type: String,
    default: 'var(--inactive-grey)'
  }
});

defineEmits(['update:modelValue']);
</script>

<style scoped>
.radio_wrap {
  display: flex;
  gap: 16px;
  flex-wrap: wrap;
}

.radio-label {
  display: flex;
  align-items: center;
  gap: 6px;
  cursor: pointer;
  user-select: none;
}
.radio-label.selected {
  background-color: rgba(0, 0, 0, 0.05); /* or use a theme color */
  padding: 8px 12px;
  border-radius: 8px;
  transition: background-color 0.3s ease;
}


.radio-input {
  display: none;
}

.radio-button {
  width: 16px;
  height: 16px;
  border-radius: 50%;
  border: 2px solid var(--inactive-grey);
  position: relative;
  transition: border-color 0.3s ease;
  display: flex;
  align-items: center;
  justify-content: center;
}

.radio-input:checked + .radio-button::after {
  content: "";
  width: 8px;
  height: 8px;
  background-color: var(--black-soft);
  border-radius: 50%;
}

.radio-text {
  color: var(--inactive-grey);
  transition: color 0.3s ease;
}

.radio-text.ractive {
  color: var(--black-soft);
}
</style>
