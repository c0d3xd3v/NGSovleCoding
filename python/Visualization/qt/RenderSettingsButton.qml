import QtQuick
import QtQuick.Controls

ToolButton {
    id: toolButton1
    text: qsTr("Tool Button")
    display: AbstractButton.TextBesideIcon
    onClicked: function() {
        console.log(swipeView.currentIndex)
        loader.source = "RenderSettings.qml"
        swipeView.currentIndex = 1
    }
}
