<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="3.18.3-Zürich" minScale="1e+08" styleCategories="AllStyleCategories" hasScaleBasedVisibilityFlag="0" maxScale="0">
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
    <property key="WMSBackgroundLayer" value="false"/>
    <property key="WMSPublishDataSourceUrl" value="false"/>
    <property key="embeddedWidgets/count" value="0"/>
    <property key="identify/format" value="Value"/>
  </customproperties>
  <pipe>
    <provider>
      <resampling zoomedInResamplingMethod="nearestNeighbour" maxOversampling="2" zoomedOutResamplingMethod="nearestNeighbour" enabled="false"/>
    </provider>
    <rasterrenderer classificationMax="24" band="5" type="singlebandpseudocolor" nodataColor="" opacity="1" classificationMin="-24" alphaBand="-1">
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
        <colorrampshader minimumValue="-24" classificationMode="2" labelPrecision="4" clip="0" maximumValue="24" colorRampType="DISCRETE">
          <colorramp type="gradient" name="[source]">
            <Option type="Map">
              <Option type="QString" name="color1" value="123,50,148,255"/>
              <Option type="QString" name="color2" value="0,136,55,255"/>
              <Option type="QString" name="discrete" value="0"/>
              <Option type="QString" name="rampType" value="gradient"/>
              <Option type="QString" name="stops" value="0.25;194,165,207,255:0.5;247,247,247,255:0.75;166,219,160,255"/>
            </Option>
            <prop k="color1" v="123,50,148,255"/>
            <prop k="color2" v="0,136,55,255"/>
            <prop k="discrete" v="0"/>
            <prop k="rampType" v="gradient"/>
            <prop k="stops" v="0.25;194,165,207,255:0.5;247,247,247,255:0.75;166,219,160,255"/>
          </colorramp>
          <item label="&lt;= -22,0000" color="#7b3294" value="-22" alpha="255"/>
          <item label="-22,0000 - -20,0000" color="#87469e" value="-20" alpha="255"/>
          <item label="-20,0000 - -18,0000" color="#945aa9" value="-18" alpha="255"/>
          <item label="-18,0000 - -16,0000" color="#a06eb3" value="-16" alpha="255"/>
          <item label="-16,0000 - -14,0000" color="#ad82bd" value="-14" alpha="255"/>
          <item label="-14,0000 - -12,0000" color="#b996c8" value="-12" alpha="255"/>
          <item label="-12,0000 - -10,0000" color="#c5a9d1" value="-10" alpha="255"/>
          <item label="-10,0000 - -8,0000" color="#ceb7d8" value="-8" alpha="255"/>
          <item label="-8,0000 - -6,0000" color="#d7c5df" value="-6" alpha="255"/>
          <item label="-6,0000 - -4,0000" color="#e0d4e6" value="-4" alpha="255"/>
          <item label="-4,0000 - -2,0000" color="#eae2ed" value="-2" alpha="255"/>
          <item label="-2,0000 - 0,0000" color="#f3f0f4" value="0" alpha="255"/>
          <item label="0,0000 - 2,0000" color="#f0f5f0" value="2" alpha="255"/>
          <item label="2,0000 - 4,0000" color="#e2f0e1" value="4" alpha="255"/>
          <item label="4,0000 - 6,0000" color="#d4ebd1" value="6" alpha="255"/>
          <item label="6,0000 - 8,0000" color="#c6e6c2" value="8" alpha="255"/>
          <item label="8,0000 - 10,0000" color="#b8e1b3" value="10" alpha="255"/>
          <item label="10,0000 - 12,0000" color="#aadda4" value="12" alpha="255"/>
          <item label="12,0000 - 14,0000" color="#90d092" value="14" alpha="255"/>
          <item label="14,0000 - 16,0000" color="#73c280" value="16" alpha="255"/>
          <item label="16,0000 - 18,0000" color="#56b46e" value="18" alpha="255"/>
          <item label="18,0000 - 20,0000" color="#39a55b" value="20" alpha="255"/>
          <item label="20,0000 - 22,0000" color="#1c9749" value="22" alpha="255"/>
          <item label="> 22,0000" color="#008837" value="inf" alpha="255"/>
          <rampLegendSettings maximumLabel="" suffix="" direction="0" prefix="" minimumLabel="" orientation="2" useContinuousLegend="1">
            <numericFormat id="basic">
              <Option type="Map">
                <Option type="QChar" name="decimal_separator" value=""/>
                <Option type="int" name="decimals" value="6"/>
                <Option type="int" name="rounding_type" value="0"/>
                <Option type="bool" name="show_plus" value="false"/>
                <Option type="bool" name="show_thousand_separator" value="true"/>
                <Option type="bool" name="show_trailing_zeros" value="false"/>
                <Option type="QChar" name="thousand_separator" value=""/>
              </Option>
            </numericFormat>
          </rampLegendSettings>
        </colorrampshader>
      </rastershader>
    </rasterrenderer>
    <brightnesscontrast gamma="1" brightness="0" contrast="0"/>
    <huesaturation colorizeGreen="128" colorizeRed="255" colorizeBlue="128" grayscaleMode="0" saturation="0" colorizeOn="0" colorizeStrength="100"/>
    <rasterresampler maxOversampling="2"/>
    <resamplingStage>resamplingFilter</resamplingStage>
  </pipe>
  <blendMode>0</blendMode>
</qgis>
