<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="3.18.3-Zürich" styleCategories="AllStyleCategories" maxScale="0" hasScaleBasedVisibilityFlag="0" minScale="1e+08">
  <flags>
    <Identifiable>1</Identifiable>
    <Removable>1</Removable>
    <Searchable>1</Searchable>
    <Private>0</Private>
  </flags>
  <temporal enabled="0" mode="0" fetchMode="0">
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
      <resampling zoomedOutResamplingMethod="nearestNeighbour" maxOversampling="2" enabled="false" zoomedInResamplingMethod="nearestNeighbour"/>
    </provider>
    <rasterrenderer band="-1" nodataColor="" type="singlebandpseudocolor" alphaBand="-1" classificationMax="24" classificationMin="-24" opacity="1">
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
        <colorrampshader maximumValue="24" minimumValue="-24" classificationMode="2" clip="0" labelPrecision="4" colorRampType="DISCRETE">
          <colorramp name="[source]" type="gradient">
            <Option type="Map">
              <Option name="color1" type="QString" value="123,50,148,255"/>
              <Option name="color2" type="QString" value="0,136,55,255"/>
              <Option name="discrete" type="QString" value="0"/>
              <Option name="rampType" type="QString" value="gradient"/>
              <Option name="stops" type="QString" value="0.25;194,165,207,255:0.5;247,247,247,255:0.75;166,219,160,255"/>
            </Option>
            <prop v="123,50,148,255" k="color1"/>
            <prop v="0,136,55,255" k="color2"/>
            <prop v="0" k="discrete"/>
            <prop v="gradient" k="rampType"/>
            <prop v="0.25;194,165,207,255:0.5;247,247,247,255:0.75;166,219,160,255" k="stops"/>
          </colorramp>
          <item label="&lt;= -22,0000" alpha="255" color="#7b3294" value="-22"/>
          <item label="-22,0000 - -20,0000" alpha="255" color="#87469e" value="-20"/>
          <item label="-20,0000 - -18,0000" alpha="255" color="#945aa9" value="-18"/>
          <item label="-18,0000 - -16,0000" alpha="255" color="#a06eb3" value="-16"/>
          <item label="-16,0000 - -14,0000" alpha="255" color="#ad82bd" value="-14"/>
          <item label="-14,0000 - -12,0000" alpha="255" color="#b996c8" value="-12"/>
          <item label="-12,0000 - -10,0000" alpha="255" color="#c5a9d1" value="-10"/>
          <item label="-10,0000 - -8,0000" alpha="255" color="#ceb7d8" value="-8"/>
          <item label="-8,0000 - -6,0000" alpha="255" color="#d7c5df" value="-6"/>
          <item label="-6,0000 - -4,0000" alpha="255" color="#e0d4e6" value="-4"/>
          <item label="-4,0000 - -2,0000" alpha="255" color="#eae2ed" value="-2"/>
          <item label="-2,0000 - 0,0000" alpha="255" color="#f3f0f4" value="0"/>
          <item label="0,0000 - 2,0000" alpha="255" color="#f0f5f0" value="2"/>
          <item label="2,0000 - 4,0000" alpha="255" color="#e2f0e1" value="4"/>
          <item label="4,0000 - 6,0000" alpha="255" color="#d4ebd1" value="6"/>
          <item label="6,0000 - 8,0000" alpha="255" color="#c6e6c2" value="8"/>
          <item label="8,0000 - 10,0000" alpha="255" color="#b8e1b3" value="10"/>
          <item label="10,0000 - 12,0000" alpha="255" color="#aadda4" value="12"/>
          <item label="12,0000 - 14,0000" alpha="255" color="#90d092" value="14"/>
          <item label="14,0000 - 16,0000" alpha="255" color="#73c280" value="16"/>
          <item label="16,0000 - 18,0000" alpha="255" color="#56b46e" value="18"/>
          <item label="18,0000 - 20,0000" alpha="255" color="#39a55b" value="20"/>
          <item label="20,0000 - 22,0000" alpha="255" color="#1c9749" value="22"/>
          <item label="> 22,0000" alpha="255" color="#008837" value="inf"/>
          <rampLegendSettings suffix="" useContinuousLegend="1" orientation="2" direction="0" minimumLabel="" prefix="" maximumLabel="">
            <numericFormat id="basic">
              <Option type="Map">
                <Option name="decimal_separator" type="QChar" value=""/>
                <Option name="decimals" type="int" value="6"/>
                <Option name="rounding_type" type="int" value="0"/>
                <Option name="show_plus" type="bool" value="false"/>
                <Option name="show_thousand_separator" type="bool" value="true"/>
                <Option name="show_trailing_zeros" type="bool" value="false"/>
                <Option name="thousand_separator" type="QChar" value=""/>
              </Option>
            </numericFormat>
          </rampLegendSettings>
        </colorrampshader>
      </rastershader>
    </rasterrenderer>
    <brightnesscontrast brightness="0" gamma="1" contrast="0"/>
    <huesaturation colorizeGreen="128" colorizeRed="255" saturation="0" grayscaleMode="0" colorizeStrength="100" colorizeBlue="128" colorizeOn="0"/>
    <rasterresampler maxOversampling="2"/>
    <resamplingStage>resamplingFilter</resamplingStage>
  </pipe>
  <blendMode>0</blendMode>
</qgis>
