import QtQuick
import QtQuick.Controls

ToolButton {
    id: exportSettingsButton
    text: qsTr("Tool Button")
    display: AbstractButton.TextBesideIcon
    onClicked: function() {
        loader.source = "ExportSettings.qml"
        swipeView.currentIndex = 1
    }
}
