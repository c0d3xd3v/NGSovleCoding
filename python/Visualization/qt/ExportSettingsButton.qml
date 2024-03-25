import QtQuick
import QtQuick.Controls

ToolButton {
    id: exportSettingsButton
    text: qsTr("Tool Button")
    display: AbstractButton.TextBesideIcon
    onClicked: function() {
        console.log(swipeView.currentIndex)
        loader.source = "ExportSettings.qml"
        swipeView.currentIndex = 1
    }
}
