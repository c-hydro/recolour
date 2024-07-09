<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="3.18.3-Zürich" styleCategories="AllStyleCategories" maxScale="0" hasScaleBasedVisibilityFlag="0" minScale="1e+08">
  <flags>
    <Identifiable>1</Identifiable>
    <Removable>1</Removable>
    <Searchable>1</Searchable>
    <Private>0</Private>
  </flags>
  <temporal fetchMode="0" enabled="0" mode="0">
    <fixedRange>
      <start></start>
      <end></end>
    </fixedRange>
  </temporal>
  <customproperties>
    <property value="false" key="WMSBackgroundLayer"/>
    <property value="false" key="WMSPublishDataSourceUrl"/>
    <property value="0" key="embeddedWidgets/count"/>
    <property value="Value" key="identify/format"/>
  </customproperties>
  <pipe>
    <provider>
      <resampling enabled="false" maxOversampling="2" zoomedInResamplingMethod="nearestNeighbour" zoomedOutResamplingMethod="nearestNeighbour"/>
    </provider>
    <rasterrenderer alphaBand="-1" nodataColor="" classificationMax="48" classificationMin="0" opacity="1" band="4" type="singlebandpseudocolor">
      <rasterTransparency/>
      <minMaxOrigin>
        <limits>None</limits>
        <extent>WholeRaster</extent>
        <statAccuracy>Estimated</statAccuracy>
        <cumulativeCutLower>0.02</cumulativeCutLower>
        <cumulativeCutUpper>0.98</cumulativeCutUpper>
        <stdDevFactor>2</stdDevFactor>
      </minMaxOrigin>
      <rastershader>
        <colorrampshader clip="0" maximumValue="48" colorRampType="DISCRETE" labelPrecision="0" minimumValue="0" classificationMode="2">
          <colorramp name="[source]" type="gradient">
            <Option type="Map">
              <Option name="color1" value="241,238,246,255" type="QString"/>
              <Option name="color2" value="152,0,67,255" type="QString"/>
              <Option name="discrete" value="0" type="QString"/>
              <Option name="rampType" value="gradient" type="QString"/>
              <Option name="stops" value="0.25;215,181,216,255:0.5;223,101,176,255:0.75;221,28,119,255" type="QString"/>
            </Option>
            <prop k="color1" v="241,238,246,255"/>
            <prop k="color2" v="152,0,67,255"/>
            <prop k="discrete" v="0"/>
            <prop k="rampType" v="gradient"/>
            <prop k="stops" v="0.25;215,181,216,255:0.5;223,101,176,255:0.75;221,28,119,255"/>
          </colorramp>
          <item color="#f1eef6" label="&lt;= 4" value="4" alpha="255"/>
          <item color="#e8daec" label="4 - 8" value="8" alpha="255"/>
          <item color="#dec5e1" label="8 - 12" value="12" alpha="255"/>
          <item color="#d8aed5" label="12 - 16" value="16" alpha="255"/>
          <item color="#db91c6" label="16 - 20" value="20" alpha="255"/>
          <item color="#de73b7" label="20 - 24" value="24" alpha="255"/>
          <item color="#df58a6" label="24 - 28" value="28" alpha="255"/>
          <item color="#de3d91" label="28 - 32" value="32" alpha="255"/>
          <item color="#de227c" label="32 - 36" value="36" alpha="255"/>
          <item color="#ca1469" label="36 - 40" value="40" alpha="255"/>
          <item color="#b10a56" label="40 - 44" value="44" alpha="255"/>
          <item color="#980043" label="> 44" value="inf" alpha="255"/>
          <rampLegendSettings orientation="2" direction="0" prefix="" minimumLabel="" maximumLabel="" suffix="" useContinuousLegend="1">
            <numericFormat id="basic">
              <Option type="Map">
                <Option name="decimal_separator" value="" type="QChar"/>
                <Option name="decimals" value="6" type="int"/>
                <Option name="rounding_type" value="0" type="int"/>
                <Option name="show_plus" value="false" type="bool"/>
                <Option name="show_thousand_separator" value="true" type="bool"/>
                <Option name="show_trailing_zeros" value="false" type="bool"/>
                <Option name="thousand_separator" value="" type="QChar"/>
              </Option>
            </numericFormat>
          </rampLegendSettings>
        </colorrampshader>
      </rastershader>
    </rasterrenderer>
    <brightnesscontrast brightness="0" contrast="0" gamma="1"/>
    <huesaturation colorizeGreen="128" colorizeOn="0" colorizeRed="255" colorizeBlue="128" saturation="0" colorizeStrength="100" grayscaleMode="0"/>
    <rasterresampler maxOversampling="2"/>
    <resamplingStage>resamplingFilter</resamplingStage>
  </pipe>
  <blendMode>0</blendMode>
</qgis>
